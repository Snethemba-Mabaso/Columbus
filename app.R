# =============================================================================
#  Columbus Spatial Statistics Dashboard
#  -----------------------------------------------------------------------------
#  Back-end data : Columbus, Ohio neighbourhood crime dataset (49 polygons)
#                  from the {spData} package (fallback: {spdep}).
#  Methods       : 1. Global Moran's I   (global spatial autocorrelation)
#                  2. Local Moran's I    (LISA clusters & outliers)
#                  3. Getis-Ord Gi*      (hot / cold spot analysis)
#
#  Install once:
#    install.packages(c("shiny","bslib","sf","spdep","spData",
#                       "ggplot2","DT","dplyr"))
#  Run:
#    shiny::runApp("app.R")
# =============================================================================

pkgs <- c("shiny", "bslib", "sf", "spdep", "spData", "ggplot2", "DT", "dplyr")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Please install missing packages: install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}

library(shiny)
library(bslib)
library(sf)
library(spdep)
library(ggplot2)
library(DT)
library(dplyr)

# spdep also exports a function called card(); make sure UI cards use bslib's
card <- bslib::card

# -----------------------------------------------------------------------------
# 1. DATA (spatial back-end)
# -----------------------------------------------------------------------------
load_columbus <- function() {
  # Search both packages for the Columbus file, whatever the format/folder
  # (location differs between package versions: .gpkg in newer, .shp in older)
  pkg_dirs <- c(system.file(package = "spData"), system.file(package = "spdep"))
  pkg_dirs <- pkg_dirs[nzchar(pkg_dirs)]
  paths <- unlist(lapply(pkg_dirs, function(d)
    list.files(d, pattern = "^columbus\\.(gpkg|shp)$", recursive = TRUE,
               full.names = TRUE, ignore.case = TRUE)))
  if (!length(paths)) {
    stop("Could not find the Columbus spatial file. Update the packages with ",
         "install.packages(c('spData','spdep')) and try again.")
  }
  # prefer .gpkg (self-contained) over .shp
  paths <- paths[order(!grepl("\\.gpkg$", paths, ignore.case = TRUE))]
  x <- st_read(paths[1], quiet = TRUE)
  # Columbus uses arbitrary planar coordinates (no real-world CRS). Drop the
  # placeholder "engineering" CRS so ggplot2 never tries to reproject to lon/lat.
  suppressWarnings(sf::st_crs(x) <- NA)
  x <- st_make_valid(x)
  x$AREA_ID <- seq_len(nrow(x))
  x
}

columbus <- load_columbus()
coords   <- unname(suppressWarnings(st_coordinates(st_centroid(st_geometry(columbus)))))

all_vars <- c(
  "Crime rate (CRIME)"                  = "CRIME",
  "Housing value, $1000s (HOVAL)"       = "HOVAL",
  "Household income, $1000s (INC)"      = "INC",
  "Open space (OPEN)"                   = "OPEN",
  "Plumbing deficiency % (PLUMB)"       = "PLUMB",
  "Distance to CBD (DISCBD)"            = "DISCBD"
)
var_choices <- all_vars[all_vars %in% names(columbus)]
var_label   <- function(v) names(var_choices)[var_choices == v]

# -----------------------------------------------------------------------------
# 2. HELPERS
# -----------------------------------------------------------------------------
build_nb <- function(type, k) {
  switch(type,
         queen = poly2nb(columbus, queen = TRUE),
         rook  = poly2nb(columbus, queen = FALSE),
         knn   = knn2nb(knearneigh(coords, k = k)))
}

nb_edges <- function(nb) {
  do.call(rbind, lapply(seq_along(nb), function(i) {
    j <- nb[[i]]; j <- j[j > 0]
    if (!length(j)) return(NULL)
    data.frame(x = coords[i, 1], y = coords[i, 2],
               xend = coords[j, 1], yend = coords[j, 2])
  }))
}

quadrant_df <- function(x, lw) {
  z   <- as.numeric(scale(x))
  lag <- lag.listw(lw, z, zero.policy = TRUE)
  data.frame(z = z, lag = lag,
             quad = ifelse(z >= 0 & lag >= 0, "HH",
                    ifelse(z <  0 & lag <  0, "LL",
                    ifelse(z >= 0 & lag <  0, "HL", "LH"))))
}

lisa_cols <- c("High-High" = "#d7191c", "Low-Low" = "#2c7bb6",
               "High-Low" = "#fdae61", "Low-High" = "#74add1",
               "Not significant" = "#e0e0e0")
gi_cols   <- c("Hot spot" = "#d7191c", "Cold spot" = "#2c7bb6",
               "Not significant" = "#e0e0e0")

theme_map <- function() {
  theme_void(base_size = 13) +
    theme(plot.title    = element_text(face = "bold", size = 14, margin = margin(b = 4)),
          plot.subtitle = element_text(colour = "grey35", size = 10.5, margin = margin(b = 8)),
          legend.position = "bottom",
          legend.title = element_text(size = 10),
          plot.margin = margin(6, 6, 6, 6))
}
theme_chart <- function() {
  theme_minimal(base_size = 13) +
    theme(plot.title    = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(colour = "grey35", size = 10.5),
          panel.grid.minor = element_blank())
}

fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", formatC(p, digits = 3, format = "f"))

# -----------------------------------------------------------------------------
# 3. UI
# -----------------------------------------------------------------------------
theme <- bs_theme(
  version = 5, bootswatch = "flatly",
  primary = "#2c3e50", secondary = "#95a5a6",
  base_font = font_collection("system-ui", "-apple-system", "Segoe UI", "Roboto", "sans-serif"),
  "card-border-radius" = "12px"
)

ui <- page_sidebar(
  title = tags$span(icon("map-location-dot"), " Columbus Spatial Statistics Explorer"),
  theme = theme,
  fillable = FALSE,
  tags$head(tags$style(HTML("
    .bslib-value-box .value-box-value { font-size: 1.6rem; }
    .sidebar-section { font-size: .75rem; letter-spacing: .08em; text-transform: uppercase;
                       color: #7b8a8b; margin: 6px 0 4px; font-weight: 600; }
    .interp { background:#f4f7f9; border-left:4px solid #2c3e50; padding:12px 16px;
              border-radius:6px; font-size:.95rem; }
  "))),
  sidebar = sidebar(
    width = 330,
    div(class = "sidebar-section", "1 · Data"),
    selectInput("var", "Variable to analyse", choices = var_choices, selected = "CRIME"),

    div(class = "sidebar-section", "2 · Spatial weights"),
    selectInput("nb_type", "Neighbour definition",
                choices = c("Queen contiguity" = "queen",
                            "Rook contiguity"  = "rook",
                            "K nearest neighbours" = "knn")),
    conditionalPanel("input.nb_type == 'knn'",
                     sliderInput("k", "Number of neighbours (k)", 2, 10, 4, step = 1)),
    selectInput("style", "Weights standardisation",
                choices = c("Row-standardised (W)" = "W", "Binary (B)" = "B")),

    div(class = "sidebar-section", "3 · Method"),
    radioButtons("method", NULL,
                 choices = c("Global Moran's I"          = "global",
                             "Local Moran's I (LISA)"    = "lisa",
                             "Getis-Ord Gi* hot spots"   = "gi")),
    sliderInput("alpha", "Significance level (α)", 0.01, 0.10, 0.05, step = 0.01),
    conditionalPanel("input.method == 'global'",
                     sliderInput("nsim", "Monte Carlo simulations", 99, 9999, 999, step = 100)),
    conditionalPanel("input.method != 'global'",
                     selectInput("padj", "Multiple-testing correction",
                                 choices = c("None" = "none", "Bonferroni" = "bonferroni",
                                             "Benjamini-Hochberg (FDR)" = "BH"))),
    actionButton("run", "Calculate", icon = icon("play"),
                 class = "btn-primary w-100 mt-2"),
    tags$small(class = "text-muted d-block mt-2",
               "Change any setting, then press Calculate to refresh results.")
  ),

  navset_card_tab(
    id = "tabs",
    nav_panel("Explore", icon = icon("compass"),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("Spatial distribution"), plotOutput("map_var", height = 420)),
        card(card_header("Neighbourhood structure"), plotOutput("map_nb", height = 420))
      ),
      layout_columns(
        col_widths = c(7, 5),
        card(card_header("Distribution"), plotOutput("hist_var", height = 260)),
        card(card_header("Summary statistics"), tableOutput("summary_tbl"))
      )
    ),
    nav_panel("Analysis", icon = icon("chart-line"),
      uiOutput("value_boxes"),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header(textOutput("plot1_title")), plotOutput("plot1", height = 430)),
        card(card_header(textOutput("plot2_title")), plotOutput("plot2", height = 430))
      ),
      card(card_header("Interpretation"), uiOutput("interpretation")),
      card(card_header(div(class = "d-flex justify-content-between align-items-center",
                           "Detailed results",
                           downloadButton("dl", "Download CSV", class = "btn-sm btn-outline-primary"))),
           DTOutput("results_tbl"))
    ),
    nav_panel("Data", icon = icon("table"),
      card(DTOutput("data_tbl"))
    ),
    nav_panel("About", icon = icon("circle-info"),
      card(
        h4("About this dashboard"),
        p("This app explores the classic ", strong("Columbus, Ohio"),
          " dataset (Anselin, 1988): 49 neighbourhoods with crime, housing value, income and other indicators."),
        h5("Methods"),
        tags$ul(
          tags$li(strong("Global Moran's I: "), "tests whether similar values cluster across the whole study area. ",
                  "I > 0 indicates positive spatial autocorrelation, I < 0 dispersion. Both an analytical (normal approximation) ",
                  "and a Monte Carlo permutation p-value are reported."),
          tags$li(strong("Local Moran's I (LISA): "), "decomposes the global statistic into local contributions, identifying ",
                  "High-High and Low-Low clusters and High-Low / Low-High spatial outliers."),
          tags$li(strong("Getis-Ord Gi*: "), "compares each neighbourhood and its neighbours with the global mean, ",
                  "flagging statistically significant hot spots (high values) and cold spots (low values).")
        ),
        h5("Tips"),
        tags$ul(
          tags$li("Row-standardised weights make Moran's I comparable across neighbour definitions."),
          tags$li("With 49 areas and many local tests, consider a multiple-testing correction (Bonferroni is conservative, FDR less so)."),
          tags$li("Results are sensitive to the choice of spatial weights; compare several definitions.")
        )
      )
    )
  )
)

# -----------------------------------------------------------------------------
# 4. SERVER
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  nb <- reactive(build_nb(input$nb_type, input$k))
  lw <- reactive(nb2listw(nb(), style = input$style, zero.policy = TRUE))

  # ---- Explore tab ----------------------------------------------------------
  output$map_var <- renderPlot({
    ggplot(columbus) +
      geom_sf(aes(fill = .data[[input$var]]), colour = "white", linewidth = 0.3) +
      scale_fill_viridis_c(option = "magma", direction = -1, name = var_label(input$var)) +
      labs(title = var_label(input$var), subtitle = "Columbus, Ohio neighbourhoods") +
      guides(fill = guide_colourbar(barwidth = 14, barheight = 0.6, title.position = "top")) +
      coord_sf(datum = NA) +
      theme_map()
  })

  output$map_nb <- renderPlot({
    nbl <- nb(); e <- nb_edges(nbl)
    ggplot() +
      geom_sf(data = columbus, fill = "#f4f7f9", colour = "grey70", linewidth = 0.3) +
      geom_segment(data = e, aes(x = x, y = y, xend = xend, yend = yend),
                   colour = "#2c3e50", alpha = 0.6, linewidth = 0.4) +
      geom_point(data = data.frame(x = coords[, 1], y = coords[, 2]),
                 aes(x, y), colour = "#d7191c", size = 1.6) +
      labs(title = "Neighbour links",
           subtitle = sprintf("Mean links per area: %.2f  |  min: %d  |  max: %d",
                              mean(spdep::card(nbl)), min(spdep::card(nbl)), max(spdep::card(nbl)))) +
      coord_sf(datum = NA) +
      theme_map()
  })

  output$hist_var <- renderPlot({
    x <- columbus[[input$var]]
    ggplot(data.frame(x = x), aes(x)) +
      geom_histogram(bins = 15, fill = "#2c3e50", colour = "white", alpha = 0.9) +
      geom_vline(xintercept = mean(x), colour = "#d7191c", linewidth = 0.9, linetype = "dashed") +
      labs(x = var_label(input$var), y = "Count",
           subtitle = "Dashed line = mean") +
      theme_chart()
  })

  output$summary_tbl <- renderTable({
    x <- columbus[[input$var]]
    data.frame(Statistic = c("N", "Mean", "Median", "Std. dev.", "Min", "Max"),
               Value = c(length(x), mean(x), median(x), sd(x), min(x), max(x)) |>
                 round(2) |> format(nsmall = 2))
  }, striped = TRUE, width = "100%")

  # ---- Analysis (runs on button press) --------------------------------------
  analysis <- eventReactive(input$run, {
    x      <- columbus[[input$var]]
    w      <- lw()
    method <- input$method
    alpha  <- input$alpha
    out <- list(method = method, var = input$var, alpha = alpha,
                nb_type = input$nb_type, style = input$style)

    if (method == "global") {
      mt <- moran.test(x, w, zero.policy = TRUE, alternative = "two.sided")
      mc <- moran.mc(x, w, nsim = input$nsim, zero.policy = TRUE, alternative = "two.sided")
      out$I    <- unname(mt$estimate["Moran I statistic"])
      out$EI   <- unname(mt$estimate["Expectation"])
      out$VarI <- unname(mt$estimate["Variance"])
      out$z    <- unname(mt$statistic)
      out$p    <- mt$p.value
      out$p_mc <- mc$p.value
      out$mc   <- mc$res
      out$nsim <- input$nsim
      out$scatter <- quadrant_df(x, w)
      out$table <- data.frame(
        Statistic = c("Moran's I", "Expected I", "Variance", "Z-score",
                      "p-value (analytical)", "p-value (Monte Carlo)", "Simulations"),
        Value = c(round(c(out$I, out$EI, out$VarI, out$z), 4),
                  round(out$p, 4), round(out$p_mc, 4), out$nsim))
    }

    if (method == "lisa") {
      lm_res <- localmoran(x, w, zero.policy = TRUE)
      pcol   <- grep("^Pr", colnames(lm_res))[1]
      p_raw  <- lm_res[, pcol]
      p_adj  <- p.adjust(p_raw, method = input$padj)
      sc     <- quadrant_df(x, w)
      cl <- ifelse(p_adj >= alpha, "Not significant",
            ifelse(sc$quad == "HH", "High-High",
            ifelse(sc$quad == "LL", "Low-Low",
            ifelse(sc$quad == "HL", "High-Low", "Low-High"))))
      cl <- factor(cl, levels = names(lisa_cols))
      out$padj_method <- input$padj
      out$cluster <- cl
      out$scatter <- transform(sc, cluster = cl)
      out$table <- data.frame(
        Area = columbus$AREA_ID, Value = round(x, 2),
        Local_I = round(lm_res[, 1], 4), Z = round(lm_res[, 4], 3),
        p_raw = round(p_raw, 4), p_adj = round(p_adj, 4),
        Cluster = as.character(cl))
      out$Iglobal <- sum(lm_res[, 1]) / sum(unlist(w$weights)) * length(x)
    }

    if (method == "gi") {
      w_self <- nb2listw(include.self(nb()), style = input$style, zero.policy = TRUE)
      g      <- as.numeric(localG(x, w_self, zero.policy = TRUE))
      p_raw  <- 2 * pnorm(-abs(g))
      p_adj  <- p.adjust(p_raw, method = input$padj)
      cl <- ifelse(p_adj >= alpha, "Not significant",
            ifelse(g > 0, "Hot spot", "Cold spot"))
      cl <- factor(cl, levels = names(gi_cols))
      out$padj_method <- input$padj
      out$g <- g
      out$cluster <- cl
      out$table <- data.frame(
        Area = columbus$AREA_ID, Value = round(x, 2),
        Gi_star_Z = round(g, 3), p_raw = round(p_raw, 4),
        p_adj = round(p_adj, 4), Class = as.character(cl))
    }
    out
  }, ignoreNULL = FALSE)

  # ---- Value boxes -----------------------------------------------------------
  output$value_boxes <- renderUI({
    a <- analysis()
    if (a$method == "global") {
      sig <- a$p_mc < a$alpha
      layout_columns(fill = FALSE, col_widths = c(3, 3, 3, 3),
        value_box("Moran's I", sprintf("%.3f", a$I), showcase = icon("bullseye"), theme = "primary"),
        value_box("Z-score", sprintf("%.2f", a$z), showcase = icon("calculator"), theme = "secondary"),
        value_box("p-value (Monte Carlo)", fmt_p(a$p_mc), showcase = icon("dice"), theme = "secondary"),
        value_box("Conclusion", if (sig) "Significant" else "Not significant",
                  showcase = icon(if (sig) "circle-check" else "circle-xmark"),
                  theme = if (sig) "success" else "danger"))
    } else {
      cl <- a$cluster
      n_sig <- sum(cl != "Not significant")
      if (a$method == "lisa") {
        layout_columns(fill = FALSE, col_widths = c(3, 3, 3, 3),
          value_box("Significant areas", paste(n_sig, "of", length(cl)), showcase = icon("location-dot"), theme = "primary"),
          value_box("High-High clusters", sum(cl == "High-High"), showcase = icon("fire"), theme = "danger"),
          value_box("Low-Low clusters", sum(cl == "Low-Low"), showcase = icon("snowflake"), theme = "info"),
          value_box("Spatial outliers", sum(cl %in% c("High-Low", "Low-High")), showcase = icon("shuffle"), theme = "warning"))
      } else {
        layout_columns(fill = FALSE, col_widths = c(4, 4, 4),
          value_box("Significant areas", paste(n_sig, "of", length(cl)), showcase = icon("location-dot"), theme = "primary"),
          value_box("Hot spots", sum(cl == "Hot spot"), showcase = icon("fire"), theme = "danger"),
          value_box("Cold spots", sum(cl == "Cold spot"), showcase = icon("snowflake"), theme = "info"))
      }
    }
  })

  # ---- Plot titles -----------------------------------------------------------
  output$plot1_title <- renderText({
    switch(analysis()$method,
           global = "Moran scatterplot", lisa = "LISA cluster map", gi = "Gi* hot / cold spot map")
  })
  output$plot2_title <- renderText({
    switch(analysis()$method,
           global = "Monte Carlo reference distribution",
           lisa = "Moran scatterplot by cluster type",
           gi = "Gi* z-scores")
  })

  # ---- Plot 1 ----------------------------------------------------------------
  output$plot1 <- renderPlot({
    a <- analysis()
    if (a$method == "global") {
      sc <- a$scatter
      quad_cols <- c(HH = "#d7191c", LL = "#2c7bb6", HL = "#fdae61", LH = "#74add1")
      ggplot(sc, aes(z, lag)) +
        geom_hline(yintercept = 0, colour = "grey60") +
        geom_vline(xintercept = 0, colour = "grey60") +
        geom_point(aes(colour = quad), size = 3, alpha = 0.85) +
        geom_smooth(method = "lm", se = FALSE, colour = "#2c3e50", linewidth = 1) +
        scale_colour_manual(values = quad_cols, name = "Quadrant") +
        labs(x = paste(var_label(a$var), "(standardised)"),
             y = "Spatial lag (standardised)",
             subtitle = sprintf("Slope = Moran's I = %.3f", a$I)) +
        theme_chart() + theme(legend.position = "bottom")
    } else if (a$method == "lisa") {
      d <- columbus; d$cl <- a$cluster
      ggplot(d) +
        geom_sf(aes(fill = cl), colour = "white", linewidth = 0.3) +
        scale_fill_manual(values = lisa_cols, drop = FALSE, name = NULL) +
        labs(title = var_label(a$var),
             subtitle = sprintf("α = %.2f, correction: %s", a$alpha, a$padj_method)) +
        coord_sf(datum = NA) +
      theme_map()
    } else {
      d <- columbus; d$cl <- a$cluster
      ggplot(d) +
        geom_sf(aes(fill = cl), colour = "white", linewidth = 0.3) +
        scale_fill_manual(values = gi_cols, drop = FALSE, name = NULL) +
        labs(title = var_label(a$var),
             subtitle = sprintf("α = %.2f, correction: %s", a$alpha, a$padj_method)) +
        coord_sf(datum = NA) +
      theme_map()
    }
  })

  # ---- Plot 2 ----------------------------------------------------------------
  output$plot2 <- renderPlot({
    a <- analysis()
    if (a$method == "global") {
      ggplot(data.frame(r = a$mc), aes(r)) +
        geom_histogram(bins = 30, fill = "#bdc3c7", colour = "white") +
        geom_vline(xintercept = a$I, colour = "#d7191c", linewidth = 1.2) +
        labs(x = "Moran's I under random permutation", y = "Frequency",
             subtitle = sprintf("Observed I = %.3f (red line) | %d permutations | p = %s",
                                a$I, a$nsim, fmt_p(a$p_mc))) +
        theme_chart()
    } else if (a$method == "lisa") {
      sc <- a$scatter
      ggplot(sc, aes(z, lag, colour = cluster)) +
        geom_hline(yintercept = 0, colour = "grey60") +
        geom_vline(xintercept = 0, colour = "grey60") +
        geom_point(size = 3, alpha = 0.9) +
        scale_colour_manual(values = lisa_cols, drop = FALSE, name = NULL) +
        labs(x = paste(var_label(a$var), "(standardised)"), y = "Spatial lag (standardised)") +
        theme_chart() + theme(legend.position = "bottom")
    } else {
      d <- columbus; d$g <- a$g
      lim <- max(abs(d$g))
      ggplot(d) +
        geom_sf(aes(fill = g), colour = "white", linewidth = 0.3) +
        scale_fill_gradient2(low = "#2c7bb6", mid = "#f7f7f7", high = "#d7191c",
                             midpoint = 0, limits = c(-lim, lim), name = "Gi* z") +
        guides(fill = guide_colourbar(barwidth = 14, barheight = 0.6, title.position = "top")) +
        labs(subtitle = "|z| > 1.96 ≈ significant at 5% (before correction)") +
        coord_sf(datum = NA) +
      theme_map()
    }
  })

  # ---- Interpretation --------------------------------------------------------
  output$interpretation <- renderUI({
    a <- analysis(); v <- var_label(a$var)
    txt <- if (a$method == "global") {
      sig <- a$p_mc < a$alpha
      dir <- if (a$I > a$EI) "positive spatial autocorrelation (similar values tend to be near each other)"
             else "negative spatial autocorrelation (dissimilar values tend to be neighbours)"
      sprintf("Moran's I for <b>%s</b> is <b>%.3f</b> (expected under randomness: %.3f), z = %.2f. The Monte Carlo p-value is <b>%s</b>, so at α = %.2f we %s the null hypothesis of spatial randomness. %s",
              v, a$I, a$EI, a$z, fmt_p(a$p_mc), a$alpha,
              if (sig) "<b>reject</b>" else "<b>fail to reject</b>",
              if (sig) paste("This indicates", dir, ".") else "There is no strong evidence of global spatial structure.")
    } else if (a$method == "lisa") {
      cl <- a$cluster
      sprintf("Of 49 neighbourhoods, <b>%d</b> show significant local spatial association for <b>%s</b> (α = %.2f): %d High-High, %d Low-Low, %d High-Low and %d Low-High. High-High / Low-Low areas are clusters of similar values; High-Low / Low-High areas are spatial outliers that differ from their surroundings.",
              sum(cl != "Not significant"), v, a$alpha, sum(cl == "High-High"),
              sum(cl == "Low-Low"), sum(cl == "High-Low"), sum(cl == "Low-High"))
    } else {
      cl <- a$cluster
      sprintf("Getis-Ord Gi* identifies <b>%d hot spot(s)</b> and <b>%d cold spot(s)</b> for <b>%s</b> at α = %.2f (correction: %s). Hot spots are areas where high values and their neighbours' high values are jointly higher than expected by chance; cold spots are the opposite.",
              sum(cl == "Hot spot"), sum(cl == "Cold spot"), v, a$alpha, a$padj_method)
    }
    div(class = "interp", HTML(txt))
  })

  # ---- Tables ----------------------------------------------------------------
  output$results_tbl <- renderDT({
    datatable(analysis()$table, rownames = FALSE,
              options = list(pageLength = 10, dom = "tip", scrollX = TRUE))
  })

  output$dl <- downloadHandler(
    filename = function() paste0("columbus_", analysis()$method, "_", analysis()$var, ".csv"),
    content  = function(file) write.csv(analysis()$table, file, row.names = FALSE)
  )

  output$data_tbl <- renderDT({
    datatable(st_drop_geometry(columbus), rownames = FALSE,
              options = list(pageLength = 15, scrollX = TRUE))
  })
}

shinyApp(ui, server)

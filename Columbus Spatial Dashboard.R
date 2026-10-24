# =============================================================================
#  Columbus Spatial Statistics Dashboard
#  -----------------------------------------------------------------------------
#  Back-end data : Columbus, Ohio neighbourhood crime dataset (49 polygons)
#                  from the {spData} package (fallback: {spdep}).
#
#  Tabs
#    Explore    : research question, about the data, maps & descriptive stats
#    Analysis   : Global Moran's I | Local Moran's I (LISA) | Geary's C
#    Regression : OLS vs Spatial Lag (SLR) vs Spatial Error (SEM) models
#    Data       : attribute table
#
#  Install once:
#    install.packages(c("shiny","bslib","sf","spdep","spData","spatialreg",
#                       "ggplot2","DT","dplyr","leaflet"))
#  Run:
#    shiny::runApp("app.R")
# =============================================================================

pkgs <- c("shiny", "bslib", "sf", "spdep", "spData", "spatialreg",
          "ggplot2", "DT", "dplyr", "leaflet")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) {
  stop("Please install missing packages: install.packages(c(",
       paste0('"', missing_pkgs, '"', collapse = ", "), "))")
}

library(shiny)
library(bslib)
library(sf)
library(spdep)
library(spatialreg)
library(ggplot2)
library(DT)
library(dplyr)
library(leaflet)

# spdep also exports a function called card(); make sure UI cards use bslib's
card <- bslib::card

# -----------------------------------------------------------------------------
# 1. DATA (spatial back-end)
# -----------------------------------------------------------------------------
load_columbus <- function() {
  # Search both packages for the Columbus file, whatever the format/folder
  pkg_dirs <- c(system.file(package = "spData"), system.file(package = "spdep"))
  pkg_dirs <- pkg_dirs[nzchar(pkg_dirs)]
  paths <- unlist(lapply(pkg_dirs, function(d)
    list.files(d, pattern = "^columbus\\.(gpkg|shp)$", recursive = TRUE,
               full.names = TRUE, ignore.case = TRUE)))
  if (!length(paths)) {
    stop("Could not find the Columbus spatial file. Update the packages with ",
         "install.packages(c('spData','spdep')) and try again.")
  }
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
col_df   <- st_drop_geometry(columbus)

# The Columbus file has no real-world coordinates, but leaflet needs lon/lat.
# For DISPLAY ONLY on the Explore map, make a rescaled copy of the polygons.
# All statistics use the original geometry above.
georeference_columbus <- function(x, centre = c(-83.00, 39.99), km_per_unit = 2.7) {
  x  <- suppressWarnings(st_cast(x, "MULTIPOLYGON"))
  bb <- st_bbox(x)
  cx <- unname((bb["xmin"] + bb["xmax"]) / 2)
  cy <- unname((bb["ymin"] + bb["ymax"]) / 2)
  kx <- km_per_unit / (111.32 * cos(centre[2] * pi / 180))   # degrees lon per unit
  ky <- km_per_unit / 110.57                                  # degrees lat per unit
  f  <- function(m) cbind((m[, 1] - cx) * kx + centre[1],
                          (m[, 2] - cy) * ky + centre[2])
  conv <- function(geom) {
    st_multipolygon(lapply(unclass(geom), function(poly) lapply(poly, f)))
  }
  st_geometry(x) <- st_sfc(lapply(st_geometry(x), conv), crs = 4326)
  x
}
columbus_ll <- georeference_columbus(columbus)

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

# The dashboard focuses on three variables
main_vars <- var_choices[var_choices %in% c("CRIME", "HOVAL", "INC")]   # explore / analysis / dependent
x_vars    <- var_choices[var_choices %in% c("HOVAL", "INC")]            # regression explanatory variables
short_label <- function(v) trimws(sub("\\s*\\(.*$", "", sub(",.*$", "", var_label(v))))
cap <- function(x) paste0(toupper(substr(x, 1, 1)), substring(x, 2))

# -----------------------------------------------------------------------------
# 2. HELPERS
# -----------------------------------------------------------------------------
build_nb <- function(type, k) {
  switch(type,
         queen = poly2nb(columbus, queen = TRUE),
         rook  = poly2nb(columbus, queen = FALSE),
         knn   = knn2nb(knearneigh(coords, k = k)))
}

build_lw <- function(type, k, style) {
  nb2listw(build_nb(type, k), style = style, zero.policy = TRUE)
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

theme_map <- function() {
  theme_void(base_size = 13) +
    theme(plot.title    = element_text(face = "bold", size = 14, colour = "#1d2b45", margin = margin(b = 4)),
          plot.subtitle = element_text(colour = "grey35", size = 10.5, margin = margin(b = 8)),
          legend.position = "bottom",
          legend.title = element_text(size = 10),
          plot.margin = margin(6, 6, 6, 6))
}
theme_chart <- function() {
  theme_minimal(base_size = 13) +
    theme(plot.title    = element_text(face = "bold", size = 14, colour = "#1d2b45"),
          plot.subtitle = element_text(colour = "grey35", size = 10.5),
          panel.grid.minor = element_blank())
}

fmt_p <- function(p) ifelse(is.na(p), "NA",
                            ifelse(p < 0.001, "< 0.001", formatC(p, digits = 3, format = "f")))
stars <- function(p) ifelse(is.na(p), "",
                     ifelse(p < 0.001, "***",
                     ifelse(p < 0.01,  "**",
                     ifelse(p < 0.05,  "*",
                     ifelse(p < 0.1,   ".", "")))))

# Standard-normal reference plot for the z-test of a global statistic
z_plot <- function(z, alpha) {
  zc  <- qnorm(1 - alpha / 2)
  lim <- max(4, abs(z) + 0.5)
  xs  <- seq(-lim, lim, length.out = 500)
  d   <- data.frame(x = xs, y = dnorm(xs))
  ggplot(d, aes(x, y)) +
    geom_area(data = d[d$x <= -zc, ], fill = "#f5b7b1") +
    geom_area(data = d[d$x >=  zc, ], fill = "#f5b7b1") +
    geom_line(colour = "#1d2b45", linewidth = 1) +
    geom_vline(xintercept = z, colour = "#d7191c", linewidth = 1.2) +
    labs(x = "z-score", y = "Density",
         subtitle = sprintf("Observed z = %.2f (red line). Shaded = rejection region at \u03b1 = %.2f",
                            z, alpha)) +
    theme_chart()
}

# ---- Regression helpers -----------------------------------------------------
model_labels <- c(ols = "OLS (non-spatial)",
                  slm = "Spatial lag model (SLR)",
                  sem = "Spatial error model (SEM)")

# The LM test function has lived in different packages/names over versions
get_lm_tests_fn <- function() {
  for (pkg in c("spdep", "spatialreg")) {
    for (fn in c("lm.RStests", "lm.LMtests")) {
      ns <- asNamespace(pkg)
      if (exists(fn, envir = ns, inherits = FALSE)) return(get(fn, envir = ns))
    }
  }
  NULL
}

run_lm_tests <- function(ols, w) {
  fn <- get_lm_tests_fn()
  if (is.null(fn)) return(NULL)
  r <- fn(ols, w, test = c("LMerr", "LMlag", "RLMerr", "RLMlag"), zero.policy = TRUE)
  labs <- c(LMerr = "LM-Error", LMlag = "LM-Lag",
            RLMerr = "Robust LM-Error", RLMlag = "Robust LM-Lag")
  keys <- names(r)
  data.frame(
    Key       = keys,
    Test      = ifelse(keys %in% names(labs), labs[keys], keys),
    Statistic = vapply(r, function(h) unname(as.numeric(h$statistic)), numeric(1)),
    df        = vapply(r, function(h) unname(as.numeric(h$parameter)), numeric(1)),
    p         = vapply(r, function(h) as.numeric(h$p.value), numeric(1)),
    row.names = NULL, stringsAsFactors = FALSE
  )
}

decision_text <- function(t, alpha) {
  p <- setNames(t$p, t$Key)
  sig <- function(k) k %in% names(p) && !is.na(p[[k]]) && p[[k]] < alpha
  rec <- if (!sig("LMerr") && !sig("LMlag")) {
    "Neither LM test is significant, so OLS is adequate."
  } else if (sig("LMlag") && !sig("LMerr")) {
    "Only the LM-Lag test is significant, pointing to the <b>spatial lag model</b>."
  } else if (sig("LMerr") && !sig("LMlag")) {
    "Only the LM-Error test is significant, pointing to the <b>spatial error model</b>."
  } else if (sig("RLMlag") && (!sig("RLMerr") || p[["RLMlag"]] < p[["RLMerr"]])) {
    "Both LM tests are significant; the robust tests favour the <b>spatial lag model</b>."
  } else if (sig("RLMerr")) {
    "Both LM tests are significant; the robust tests favour the <b>spatial error model</b>."
  } else {
    "Both LM tests are significant but the robust versions are not decisive; compare the models by AIC."
  }
  paste("<b>Lagrange multiplier diagnostics:</b>", rec)
}

summarise_model <- function(fit, type, w) {
  if (inherits(fit, "error")) return(list(error = conditionMessage(fit)))
  tryCatch({
    if (type == "ols") {
      cf  <- coef(summary(fit))
      tab <- data.frame(Term = rownames(cf), Estimate = cf[, 1], SE = cf[, 2],
                        Stat = cf[, 3], p = cf[, 4], row.names = NULL)
      param <- summary(fit)$adj.r.squared
      param_lab <- "Adjusted R\u00b2"
    } else {
      sm  <- summary(fit)
      cf  <- sm$Coef
      tab <- data.frame(Term = rownames(cf), Estimate = cf[, 1], SE = cf[, 2],
                        Stat = cf[, 3], p = cf[, 4], row.names = NULL)
      if (type == "slm") {
        est <- fit$rho;    se <- if (is.null(fit$rho.se)) NA_real_ else fit$rho.se
        nm <- "rho (spatial lag)";     param_lab <- "\u03c1 (rho)"
      } else {
        est <- fit$lambda; se <- if (is.null(fit$lambda.se)) NA_real_ else fit$lambda.se
        nm <- "lambda (spatial error)"; param_lab <- "\u03bb (lambda)"
      }
      z <- est / se
      tab <- rbind(tab, data.frame(Term = nm, Estimate = est, SE = se,
                                   Stat = z, p = 2 * pnorm(-abs(z))))
      param <- est
    }
    res <- as.numeric(residuals(fit))
    mt  <- moran.test(res, w, zero.policy = TRUE, alternative = "two.sided")
    list(fit = fit, table = tab, aic = AIC(fit), ll = as.numeric(logLik(fit)),
         res = res, mI = unname(mt$estimate["Moran I statistic"]),
         mP = mt$p.value, param = param, param_lab = param_lab)
  }, error = function(e) list(error = conditionMessage(e)))
}

# -----------------------------------------------------------------------------
# 3. UI
# -----------------------------------------------------------------------------
theme <- bs_theme(
  version   = 5,
  bg        = "#ffffff",
  fg        = "#243047",
  primary   = "#1d2b45",   # navy
  secondary = "#6b7a8c",   # slate grey
  success   = "#2a9d8f",   # teal
  info      = "#4a7ba7",   # steel blue
  warning   = "#d99a3d",   # muted amber
  danger    = "#c75d5d",   # muted red
  base_font = font_collection("Inter", "system-ui", "-apple-system", "Segoe UI", "Roboto", "sans-serif"),
  "border-radius"      = "0.6rem",
  "card-border-radius" = "0.9rem"
)

app_css <- "
  /* ---------- page ---------- */
  body { background: #f1f4f8; }
  .container-fluid, .bslib-page-navbar > .container-fluid { padding-top: .25rem; }

  /* ---------- navbar ---------- */
  .navbar { position: sticky; top: 0; z-index: 1030;   /* frozen header */
            box-shadow: 0 2px 10px rgba(15, 25, 45, .25); padding-top: .55rem; padding-bottom: .55rem; }
  .navbar-brand { font-weight: 700; letter-spacing: .01em; color: #ffffff !important; }
  .navbar .nav-link { color: rgba(255,255,255,.82) !important; font-weight: 500;
                      padding: .5rem .9rem; border-bottom: 2px solid transparent; margin: 0 .1rem; }
  .navbar .nav-link:hover { color: #ffffff !important; }
  .navbar .nav-link.active { color: #7fd6cc !important; border-bottom-color: #7fd6cc; background: transparent !important; }
  .nav-note { color: rgba(255,255,255,.55); font-size: .82rem; padding-right: .5rem; }

  /* ---------- leaflet (no basemap) ---------- */
  .leaflet-container { background: #f4f7fb; font-family: inherit; }

  /* ---------- cards ---------- */
  .card { border: 1px solid #dfe5ee; box-shadow: 0 2px 8px rgba(20, 33, 61, .06); }
  .card-header { background: #ffffff; color: #1d2b45; font-weight: 600; font-size: .95rem;
                 border-bottom: 1px solid #e6ebf2; }
  .card-footer { background: #f7f9fc; border-top: 1px solid #e6ebf2; }

  /* ---------- sidebar ---------- */
  .bslib-sidebar-layout > .sidebar { background: #e8edf4 !important; border-right: 1px solid #d5dde8; }
  /* frozen sidebar: stays in view under the frozen navbar while the results scroll */
  @media (min-width: 768px) {
    .bslib-sidebar-layout { overflow: visible !important; }
    .bslib-sidebar-layout > .sidebar {
      position: sticky !important; top: 72px; align-self: start;
      height: auto !important; max-height: calc(100vh - 88px); overflow-y: auto;
      border: 1px solid #d5dde8; border-radius: 12px; margin: .25rem .5rem .25rem 0;
      box-shadow: 0 2px 8px rgba(20, 33, 61, .06);
    }
  }
  .sidebar-section { font-size: .72rem; letter-spacing: .09em; text-transform: uppercase;
                     color: #4a5a73; margin: 14px 0 8px; font-weight: 700;
                     border-bottom: 1px solid #cdd6e3; padding-bottom: 4px; }
  .sidebar-section:first-child { margin-top: 4px; }
  .form-label, .control-label { font-weight: 500; color: #2f3d57; font-size: .92rem; }
  .form-select, .form-control { border-color: #cdd6e3; background-color: #fbfcfe; }
  .form-select:focus, .form-control:focus { border-color: #4a7ba7; box-shadow: 0 0 0 .2rem rgba(74,123,167,.2); }

  /* ---------- sliders ---------- */
  .irs--shiny .irs-bar, .irs--shiny .irs-bar-edge { background: #1d2b45; border-color: #1d2b45; }
  .irs--shiny .irs-single, .irs--shiny .irs-from, .irs--shiny .irs-to { background: #1d2b45; }
  .irs--shiny .irs-handle { border-color: #1d2b45; }
  .irs--shiny .irs-line { background: #cdd6e3; border-color: #cdd6e3; }

  /* ---------- buttons ---------- */
  .btn { border-radius: .6rem; font-weight: 600; }
  .btn-primary { background: #1d2b45; border-color: #1d2b45; }
  .btn-primary:hover, .btn-primary:focus { background: #2e4470; border-color: #2e4470; }
  .btn-outline-primary { color: #1d2b45; border-color: #1d2b45; }
  .btn-outline-primary:hover { background: #1d2b45; color: #fff; }

  /* ---------- value boxes ---------- */
  .bslib-value-box { border: none; box-shadow: 0 3px 10px rgba(20, 33, 61, .12); }
  .bslib-value-box .value-box-value { font-size: 1.6rem; font-weight: 700; }
  .bslib-value-box .value-box-title { font-weight: 500; opacity: .9; }

  /* ---------- text callouts ---------- */
  .tab-desc { background: #ffffff; border: 1px solid #dfe5ee; border-left: 5px solid #4a7ba7;
              border-radius: 10px; padding: 11px 16px; margin-bottom: 1rem;
              font-size: .93rem; color: #34425c; box-shadow: 0 1px 4px rgba(20,33,61,.05); }
  .tab-desc .fa, .tab-desc svg { color: #4a7ba7; margin-right: 4px; }
  .interp { background: #f4f7fb; border-left: 4px solid #1d2b45; padding: 12px 16px;
            border-radius: 8px; font-size: .93rem; color: #2b3a55; }
  .rq { background: #ffffff; border: 1px solid #d5dde8; border-left: 4px solid #2a9d8f;
        padding: 10px 14px; border-radius: 8px; font-size: .95rem; margin-bottom: 8px; color: #1d2b45; }
  .best-model { background: #e7f5f3; border-left: 5px solid #2a9d8f; padding: 12px 16px;
                border-radius: 8px; font-size: .95rem; color: #1f3b3a; }
  .data-dict li { margin-bottom: 4px; font-size: .9rem; }
  .small-note { font-size: .85rem; color: #5f6e84; }

  /* ---------- tables ---------- */
  table.table thead th, table thead th { color: #1d2b45; border-bottom: 2px solid #cdd6e3 !important; }
  .table-striped > tbody > tr:nth-of-type(odd) > * { --bs-table-accent-bg: #f4f7fb; }
  table.dataTable thead th { color: #1d2b45; }
"

ui <- page_navbar(
  title = tags$span(icon("map-location-dot"), " Columbus Spatial Statistics Explorer"),
  theme = theme,
  fillable = FALSE,
  bg = "#1d2b45",
  inverse = TRUE,
  header = tagList(
    tags$link(rel = "stylesheet",
              href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap"),
    tags$style(HTML(app_css))
  ),

  # ===========================================================================
  # EXPLORE
  # ===========================================================================
  nav_panel(
    "Explore", icon = icon("compass"),
    layout_sidebar(
      sidebar = sidebar(
        width = 400, open = "always",
        div(class = "sidebar-section", "Research question"),
        div(class = "rq",
            tags$p(strong("Is neighbourhood crime in Columbus, Ohio spatially clustered, and do housing value and household income explain it once spatial dependence is taken into account?"))),
        tags$p(class = "small-note", "This dashboard works through three connected questions:"),
        tags$ol(
          tags$li(strong("Global pattern: "), "is crime spatially autocorrelated across the city? (Moran's I, Geary's C)"),
          tags$li(strong("Local pattern: "), "where are the clusters and spatial outliers? (LISA)"),
          tags$li(strong("Explanation: "), "does a spatial regression (lag or error) describe crime better than ordinary least squares?")
        ),
        tags$p(class = "small-note",
               HTML("<b>Null hypothesis:</b> values are distributed randomly in space (no spatial autocorrelation).")),

        div(class = "sidebar-section", "About the data"),
        tags$p("The classic ", strong("Columbus, Ohio"),
               " dataset (Anselin, 1988): ", strong("49 neighbourhoods"),
               " (polygons) with crime and socio-economic indicators. Loaded from the ",
               tags$code("spData"), " package."),
        tags$ul(class = "data-dict",
          tags$li(strong("CRIME: "), "residential burglaries and vehicle thefts per 1,000 households"),
          tags$li(strong("HOVAL: "), "housing value ($1,000s)"),
          tags$li(strong("INC: "), "household income ($1,000s)"),
          tags$li(strong("OPEN: "), "amount of open space in the neighbourhood"),
          tags$li(strong("PLUMB: "), "% of housing units without plumbing"),
          tags$li(strong("DISCBD: "), "distance to the central business district")
        ),
        tags$p(class = "small-note",
               "Coordinates are planar and arbitrary (no real-world map projection), so maps show relative shape and position only.")
      ),

      # ---- main area -------------------------------------------------------
      div(class = "mb-2",
          selectInput("ex_var", "Variable to explore", choices = main_vars,
                      selected = "CRIME", width = "340px")),
      layout_columns(
        col_widths = c(7, 5),
        card(card_header("Spatial distribution (click a neighbourhood for details)"),
             leafletOutput("map_var", height = 480),
             card_footer(uiOutput("map_interp"))),
        layout_columns(
          col_widths = c(12, 12),
          card(card_header("Distribution"), plotOutput("hist_var", height = 220),
               card_footer(uiOutput("hist_interp"))),
          card(card_header("Summary statistics"), tableOutput("summary_tbl"))
        )
      ),
      card(card_header(div(class = "d-flex justify-content-between align-items-center",
                           "Data (49 Columbus neighbourhoods)",
                           downloadButton("dl_data", "Download CSV", class = "btn-sm btn-outline-primary"))),
           DTOutput("ex_data_tbl"))
    )
  ),

  # ===========================================================================
  # ANALYSIS
  # ===========================================================================
  nav_panel(
    "Analysis", icon = icon("chart-line"),
    layout_sidebar(
      sidebar = sidebar(
        width = 340, open = "always",
        div(class = "sidebar-section", "1 \u00b7 Data"),
        selectInput("var", "Variable to analyse", choices = main_vars, selected = "CRIME"),

        div(class = "sidebar-section", "2 \u00b7 Spatial weights"),
        selectInput("nb_type", "Neighbour definition",
                    choices = c("Queen contiguity" = "queen",
                                "Rook contiguity"  = "rook",
                                "K nearest neighbours" = "knn")),
        conditionalPanel("input.nb_type == 'knn'",
                         sliderInput("k", "Number of neighbours (k)", 2, 10, 4, step = 1)),

        div(class = "sidebar-section", "3 \u00b7 Method"),
        radioButtons("method", NULL,
                     choices = c("Global Moran's I"       = "global",
                                 "Local Moran's I (LISA)" = "lisa",
                                 "Geary's C"              = "geary")),
        sliderInput("alpha", "Significance level (\u03b1)", 0.01, 0.10, 0.05, step = 0.01),
        actionButton("run", "Calculate", icon = icon("play"),
                     class = "btn-primary w-100 mt-2"),
        tags$small(class = "text-muted d-block mt-2",
                   "Change any setting, then press Calculate to refresh results.")
      ),

      card(card_header("Hypothesis"), uiOutput("hypothesis")),
      uiOutput("value_boxes"),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header(textOutput("plot1_title")), plotOutput("plot1", height = 430)),
        card(card_header(textOutput("plot2_title")), plotOutput("plot2", height = 430))
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("Neighbourhood structure"), plotOutput("an_map_nb", height = 380)),
        card(card_header("Interpretation"), uiOutput("interpretation"))
      ),
      card(card_header(div(class = "d-flex justify-content-between align-items-center",
                           "Detailed results",
                           downloadButton("dl", "Download CSV", class = "btn-sm btn-outline-primary"))),
           DTOutput("results_tbl"))
    )
  ),

  # ===========================================================================
  # REGRESSION
  # ===========================================================================
  nav_panel(
    "Regression", icon = icon("square-poll-vertical"),
    layout_sidebar(
      sidebar = sidebar(
        width = 340, open = "always",
        div(class = "sidebar-section", "1 \u00b7 Model specification"),
        selectInput("reg_y", "Dependent variable", choices = main_vars, selected = "CRIME"),
        checkboxGroupInput("reg_x", "Explanatory variables", choices = x_vars,
                           selected = unname(x_vars)),

        div(class = "sidebar-section", "2 \u00b7 Spatial weights"),
        selectInput("reg_nb", "Neighbour definition",
                    choices = c("Queen contiguity" = "queen",
                                "Rook contiguity"  = "rook",
                                "K nearest neighbours" = "knn")),
        conditionalPanel("input.reg_nb == 'knn'",
                         sliderInput("reg_k", "Number of neighbours (k)", 2, 10, 4, step = 1)),

        div(class = "sidebar-section", "3 \u00b7 Model & significance"),
        radioButtons("reg_model", "Model to display",
                     choices = c("OLS (non-spatial)"        = "ols",
                                 "Spatial lag model (SLR)"   = "slm",
                                 "Spatial error model (SEM)" = "sem")),
        sliderInput("reg_alpha", "Significance level (\u03b1)", 0.01, 0.10, 0.05, step = 0.01),
        actionButton("reg_run", "Estimate models", icon = icon("play"),
                     class = "btn-primary w-100 mt-2"),
        tags$small(class = "text-muted d-block mt-2",
                   "All three models are estimated together so they can be compared. ",
                   "The radio button chooses which one is displayed.")
      ),

      uiOutput("reg_boxes"),
      card(card_header(textOutput("reg_hyp_title")), uiOutput("reg_hyp")),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header(textOutput("reg_coef_title")), DTOutput("reg_coef")),
        card(card_header("Model comparison"), tableOutput("reg_compare"),
             uiOutput("reg_cmp_best"),
             tags$p(class = "small-note mt-2",
                    "Lower AIC = better fit. Residual Moran's I tests for leftover spatial autocorrelation (p \u2265 \u03b1 is desirable)."))
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("Residual map"), plotOutput("reg_map", height = 400)),
        card(card_header("Residual Moran scatterplot"), plotOutput("reg_scatter", height = 400))
      ),
      card(card_header("Spatial dependence diagnostics (on OLS residuals)"),
           tableOutput("reg_lm")),
      card(card_header("Interpretation"), uiOutput("reg_text")),
      card(card_header("Conclusion"),
           uiOutput("reg_conclusion"))
    )
  ),

  # ===========================================================================
  # DATA
  # ===========================================================================
  nav_panel("Data", icon = icon("table"),
            card(DTOutput("data_tbl"))),

  nav_spacer(),
  nav_item(tags$span(class = "nav-note", "Spatial Statistics \u00b7 STAT7006A"))
)

# -----------------------------------------------------------------------------
# 4. SERVER
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  # ===========================================================================
  # EXPLORE
  # ===========================================================================
  output$map_var <- renderLeaflet({
    v <- input$ex_var
    x <- col_df[[v]]
    pal <- colorNumeric("magma", domain = x, reverse = TRUE)   # light = low, dark = high

    # Pop-up: all variables for the clicked neighbourhood, selected one highlighted
    popup <- vapply(seq_len(nrow(col_df)), function(i) {
      rows <- vapply(seq_along(main_vars), function(j) {
        vj <- main_vars[[j]]
        sprintf("<tr%s><td style='padding-right:14px'>%s</td><td style='text-align:right'><b>%.2f</b></td></tr>",
                if (vj == v) " style='background:#e8edf4'" else "",
                names(main_vars)[j], col_df[[vj]][i])
      }, character(1))
      sprintf("<div style='min-width:230px'><div style='font-weight:700;color:#1d2b45;margin-bottom:6px'>Neighbourhood #%d</div><table style='font-size:12px;width:100%%'>%s</table></div>",
              columbus$AREA_ID[i], paste(rows, collapse = ""))
    }, character(1))

    leaflet(columbus_ll) %>%
      addPolygons(
        fillColor = pal(x), fillOpacity = 0.8,
        color = "white", weight = 1.2,
        popup = popup, popupOptions = popupOptions(maxWidth = 320),
        label = sprintf("Neighbourhood #%d \u00b7 %s: %.2f", columbus$AREA_ID, var_label(v), x),
        highlightOptions = highlightOptions(weight = 3, color = "#1d2b45",
                                            fillOpacity = 0.95, bringToFront = TRUE)) %>%
      addLegend(position = "bottomright", pal = pal, values = x,
                title = var_label(v), opacity = 0.9)
  })

  output$hist_var <- renderPlot({
    x <- columbus[[input$ex_var]]
    ggplot(data.frame(x = x), aes(x)) +
      geom_histogram(bins = 15, fill = "#1d2b45", colour = "white", alpha = 0.9) +
      geom_vline(xintercept = mean(x), colour = "#d7191c", linewidth = 0.9, linetype = "dashed") +
      labs(x = var_label(input$ex_var), y = "Count", subtitle = "Dashed line = mean") +
      theme_chart()
  })

  output$map_interp <- renderUI({
    v <- input$ex_var; x <- columbus[[v]]; id <- columbus$AREA_ID
    f <- function(i) sprintf("#%d (%.2f)", id[i], x[i])
    top <- head(order(x, decreasing = TRUE), 3)
    bot <- head(order(x), 3)
    div(class = "interp", HTML(sprintf(
      "Darker shades mark higher values of <b>%s</b>. The highest values are in neighbourhoods %s and the lowest in %s (values range from %.2f to %.2f). Look at whether similar shades tend to sit next to each other: that visual impression is tested formally on the <b>Analysis</b> tab.",
      var_label(v), paste(vapply(top, f, character(1)), collapse = ", "),
      paste(vapply(bot, f, character(1)), collapse = ", "), min(x), max(x))))
  })

  output$hist_interp <- renderUI({
    x <- columbus[[input$ex_var]]
    m <- mean(x); md <- median(x); s <- sd(x)
    sk <- mean((x - m)^3) / (mean((x - m)^2)^1.5)
    shape <- if (sk > 0.5) "right-skewed (a few neighbourhoods have unusually high values)"
             else if (sk < -0.5) "left-skewed (a few neighbourhoods have unusually low values)"
             else "roughly symmetric"
    rel <- if (abs(m - md) < 0.1 * s) "close to" else if (m > md) "above" else "below"
    out <- boxplot.stats(x)$out
    out_txt <- if (length(out)) sprintf(
      "%d potential outlier(s) lie more than 1.5 \u00d7 IQR beyond the quartiles; outliers can strongly influence Moran's I and regression estimates.",
      length(out)) else "No outliers by the 1.5 \u00d7 IQR rule."
    div(class = "interp", HTML(sprintf(
      "The distribution of <b>%s</b> is <b>%s</b> (skewness = %.2f). The mean (%.2f) is %s the median (%.2f) and the standard deviation is %.2f. %s",
      var_label(input$ex_var), shape, sk, m, rel, md, s, out_txt)))
  })

  output$ex_data_tbl <- renderDT({
    datatable(col_df, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$dl_data <- downloadHandler(
    filename = function() "columbus_data.csv",
    content  = function(file) write.csv(col_df, file, row.names = FALSE)
  )

  output$summary_tbl <- renderTable({
    x <- columbus[[input$ex_var]]
    data.frame(Statistic = c("N", "Mean", "Median", "Std. dev.", "Min", "Max"),
               Value = format(round(c(length(x), mean(x), median(x), sd(x), min(x), max(x)), 2),
                              nsmall = 2))
  }, striped = TRUE, width = "100%")

  # ===========================================================================
  # ANALYSIS
  # ===========================================================================
  nb <- reactive(build_nb(input$nb_type, input$k))
  lw <- reactive(nb2listw(nb(), style = "W", zero.policy = TRUE))

  output$an_map_nb <- renderPlot({
    nbl <- nb(); e <- nb_edges(nbl)
    lab <- switch(input$nb_type,
                  queen = "Queen contiguity", rook = "Rook contiguity",
                  knn   = paste0(input$k, " nearest neighbours"))
    ggplot() +
      geom_sf(data = columbus, fill = "#f4f7f9", colour = "grey70", linewidth = 0.3) +
      geom_segment(data = e, aes(x = x, y = y, xend = xend, yend = yend),
                   colour = "#1d2b45", alpha = 0.6, linewidth = 0.4) +
      geom_point(data = data.frame(x = coords[, 1], y = coords[, 2]),
                 aes(x, y), colour = "#d7191c", size = 1.6) +
      labs(title = lab,
           subtitle = sprintf("Mean links per area: %.2f  |  min: %d  |  max: %d",
                              mean(spdep::card(nbl)), min(spdep::card(nbl)),
                              max(spdep::card(nbl)))) +
      coord_sf(datum = NA) +
      theme_map()
  })

  analysis <- eventReactive(input$run, {
    x      <- columbus[[input$var]]
    w      <- lw()
    method <- input$method
    alpha  <- input$alpha
    out <- list(method = method, var = input$var, alpha = alpha,
                nb_type = input$nb_type, style = "W", x = x)

    if (method == "global") {
      mt <- moran.test(x, w, zero.policy = TRUE, alternative = "two.sided")
      out$I    <- unname(mt$estimate["Moran I statistic"])
      out$EI   <- unname(mt$estimate["Expectation"])
      out$VarI <- unname(mt$estimate["Variance"])
      out$z    <- unname(mt$statistic)
      out$p    <- mt$p.value
      out$scatter <- quadrant_df(x, w)
      out$table <- data.frame(
        Statistic = c("Moran's I", "Expected I", "Variance", "Z-score", "p-value"),
        Value = c(formatC(c(out$I, out$EI, out$VarI, out$z), digits = 4, format = "f"),
                  fmt_p(out$p)))
    }

    if (method == "geary") {
      gt    <- geary.test(x, w, zero.policy = TRUE, alternative = "two.sided")
      w_row <- nb2listw(nb(), style = "W", zero.policy = TRUE)
      out$C     <- unname(gt$estimate["Geary C statistic"])
      out$EC    <- unname(gt$estimate["Expectation"])
      out$VarC  <- unname(gt$estimate["Variance"])
      out$z     <- unname(gt$statistic)
      out$p     <- gt$p.value
      out$nbmean <- lag.listw(w_row, x, zero.policy = TRUE)
      out$table <- data.frame(
        Statistic = c("Geary's C", "Expected C", "Variance", "Z-score", "p-value"),
        Value = c(formatC(c(out$C, out$EC, out$VarC, out$z), digits = 4, format = "f"),
                  fmt_p(out$p)))
    }

    if (method == "lisa") {
      lm_res <- localmoran(x, w, zero.policy = TRUE)
      pcol   <- grep("^Pr", colnames(lm_res))[1]
      p_val  <- lm_res[, pcol]
      sc     <- quadrant_df(x, w)
      cl <- ifelse(p_val >= alpha, "Not significant",
            ifelse(sc$quad == "HH", "High-High",
            ifelse(sc$quad == "LL", "Low-Low",
            ifelse(sc$quad == "HL", "High-Low", "Low-High"))))
      cl <- factor(cl, levels = names(lisa_cols))
      out$cluster <- cl
      out$scatter <- transform(sc, cluster = cl)
      out$table <- data.frame(
        Area = columbus$AREA_ID, Value = round(x, 2),
        Local_I = round(lm_res[, 1], 4), Z = round(lm_res[, 4], 3),
        p_value = round(p_val, 4), Cluster = as.character(cl))
    }
    out
  }, ignoreNULL = FALSE)

  # ---- Hypothesis -------------------------------------------------------------
  output$hypothesis <- renderUI({
    a <- analysis(); v <- var_label(a$var)
    txt <- if (a$method == "global") {
      sprintf("<b>H<sub>0</sub>:</b> <b>%s</b> is randomly distributed across Columbus neighbourhoods: there is no spatial autocorrelation (Moran's I equals its expected value of %.3f).<br><b>H<sub>1</sub>:</b> <b>%s</b> is spatially autocorrelated (Moran's I differs from %.3f).<br><b>Result:</b> p = %s, so at \u03b1 = %.2f we %s H<sub>0</sub>.",
              v, a$EI, v, a$EI, fmt_p(a$p), a$alpha,
              if (a$p < a$alpha) "<b>reject</b>" else "<b>fail to reject</b>")
    } else if (a$method == "geary") {
      sprintf("<b>H<sub>0</sub>:</b> <b>%s</b> is randomly distributed across Columbus neighbourhoods: there is no spatial autocorrelation (Geary's C equals 1).<br><b>H<sub>1</sub>:</b> <b>%s</b> is spatially autocorrelated (C differs from 1: below 1 means neighbours are similar, above 1 means they differ).<br><b>Result:</b> p = %s, so at \u03b1 = %.2f we %s H<sub>0</sub>.",
              v, v, fmt_p(a$p), a$alpha,
              if (a$p < a$alpha) "<b>reject</b>" else "<b>fail to reject</b>")
    } else {
      n_sig <- sum(a$cluster != "Not significant")
      sprintf("<b>H<sub>0</sub></b> (tested separately for each neighbourhood)<b>:</b> the value of <b>%s</b> in a neighbourhood is unrelated to the values around it (no local spatial association).<br><b>H<sub>1</sub>:</b> the neighbourhood is part of a local cluster or is a spatial outlier.<br><b>Result:</b> H<sub>0</sub> is rejected for <b>%d of 49</b> neighbourhoods at \u03b1 = %.2f (no multiple-testing correction).",
              v, n_sig, a$alpha)
    }
    div(class = "interp", HTML(txt))
  })

  # ---- Value boxes ------------------------------------------------------------
  output$value_boxes <- renderUI({
    a <- analysis()
    if (a$method %in% c("global", "geary")) {
      sig <- a$p < a$alpha
      lab <- if (a$method == "global") "Moran's I" else "Geary's C"
      val <- if (a$method == "global") a$I else a$C
      layout_columns(fill = FALSE, col_widths = c(3, 3, 3, 3),
        value_box(lab, sprintf("%.3f", val), showcase = icon("bullseye"), theme = "primary"),
        value_box("Z-score", sprintf("%.2f", a$z), showcase = icon("calculator"), theme = "secondary"),
        value_box("p-value", fmt_p(a$p), showcase = icon("flask"), theme = "secondary"),
        value_box("Conclusion", if (sig) "Significant" else "Not significant",
                  showcase = icon(if (sig) "circle-check" else "circle-xmark"),
                  theme = if (sig) "success" else "danger"))
    } else {
      cl <- a$cluster
      n_sig <- sum(cl != "Not significant")
      layout_columns(fill = FALSE, col_widths = c(3, 3, 3, 3),
        value_box("Significant areas", paste(n_sig, "of", length(cl)), showcase = icon("location-dot"), theme = "primary"),
        value_box("High-High clusters", sum(cl == "High-High"), showcase = icon("fire"), theme = "danger"),
        value_box("Low-Low clusters", sum(cl == "Low-Low"), showcase = icon("snowflake"), theme = "info"),
        value_box("Spatial outliers", sum(cl %in% c("High-Low", "Low-High")), showcase = icon("shuffle"), theme = "warning"))
    }
  })

  # ---- Plot titles ------------------------------------------------------------
  output$plot1_title <- renderText({
    switch(analysis()$method,
           global = "Moran scatterplot",
           lisa   = "LISA cluster map",
           geary  = "Z-test: observed vs. normal distribution")
  })
  output$plot2_title <- renderText({
    switch(analysis()$method,
           global = "Z-test: observed vs. normal distribution",
           lisa   = "Moran scatterplot by cluster type",
           geary  = "Value vs. neighbour average")
  })

  # ---- Plot 1 -----------------------------------------------------------------
  output$plot1 <- renderPlot({
    a <- analysis()
    if (a$method == "global") {
      sc <- a$scatter
      quad_cols <- c(HH = "#d7191c", LL = "#2c7bb6", HL = "#fdae61", LH = "#74add1")
      ggplot(sc, aes(z, lag)) +
        geom_hline(yintercept = 0, colour = "grey60") +
        geom_vline(xintercept = 0, colour = "grey60") +
        geom_point(aes(colour = quad), size = 3, alpha = 0.85) +
        geom_smooth(method = "lm", se = FALSE, colour = "#1d2b45", linewidth = 1,
                    formula = y ~ x) +
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
             subtitle = sprintf("Significance level \u03b1 = %.2f", a$alpha)) +
        coord_sf(datum = NA) +
        theme_map()
    } else {
      z_plot(a$z, a$alpha)
    }
  })

  # ---- Plot 2 -----------------------------------------------------------------
  output$plot2 <- renderPlot({
    a <- analysis()
    if (a$method == "global") {
      z_plot(a$z, a$alpha)
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
      d <- data.frame(x = a$x, nbmean = a$nbmean)
      ggplot(d, aes(x, nbmean)) +
        geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "dashed") +
        geom_point(colour = "#1d2b45", size = 3, alpha = 0.85) +
        geom_smooth(method = "lm", se = FALSE, colour = "#d7191c", linewidth = 1,
                    formula = y ~ x) +
        labs(x = var_label(a$var), y = "Average of neighbours' values",
             subtitle = "Points near the dashed 45\u00b0 line = neighbours have similar values (low C)") +
        theme_chart()
    }
  })

  # ---- Interpretation ---------------------------------------------------------
  output$interpretation <- renderUI({
    a <- analysis(); v <- var_label(a$var)
    txt <- if (a$method == "global") {
      sig <- a$p < a$alpha
      dir <- if (a$I > a$EI) "positive spatial autocorrelation (similar values tend to be near each other)"
             else "negative spatial autocorrelation (dissimilar values tend to be neighbours)"
      sprintf("Moran's I for <b>%s</b> is <b>%.3f</b> (expected under randomness: %.3f), z = %.2f. The p-value is <b>%s</b>, so at \u03b1 = %.2f we %s the null hypothesis of spatial randomness. %s",
              v, a$I, a$EI, a$z, fmt_p(a$p), a$alpha,
              if (sig) "<b>reject</b>" else "<b>fail to reject</b>",
              if (sig) paste0("This indicates ", dir, ".") else "There is no strong evidence of global spatial structure.")
    } else if (a$method == "geary") {
      sig <- a$p < a$alpha
      dir <- if (a$C < a$EC) "positive spatial autocorrelation (neighbouring values are similar)"
             else "negative spatial autocorrelation (neighbouring values differ more than expected)"
      sprintf("Geary's C for <b>%s</b> is <b>%.3f</b> (expected under randomness: %.3f), z = %.2f. The p-value is <b>%s</b>, so at \u03b1 = %.2f we %s the null hypothesis of spatial randomness. %s Note that C is inversely related to Moran's I: values below 1 mean positive autocorrelation.",
              v, a$C, a$EC, a$z, fmt_p(a$p), a$alpha,
              if (sig) "<b>reject</b>" else "<b>fail to reject</b>",
              if (sig) paste0("This indicates ", dir, ".") else "There is no strong evidence of global spatial structure.")
    } else {
      cl <- a$cluster
      sprintf("Of 49 neighbourhoods, <b>%d</b> show significant local spatial association for <b>%s</b> (\u03b1 = %.2f): %d High-High, %d Low-Low, %d High-Low and %d Low-High. High-High / Low-Low areas are clusters of similar values; High-Low / Low-High areas are spatial outliers that differ from their surroundings.",
              sum(cl != "Not significant"), v, a$alpha, sum(cl == "High-High"),
              sum(cl == "Low-Low"), sum(cl == "High-Low"), sum(cl == "Low-High"))
    }
    div(class = "interp", HTML(txt))
  })

  output$results_tbl <- renderDT({
    datatable(analysis()$table, rownames = FALSE,
              options = list(pageLength = 10, dom = "tip", scrollX = TRUE))
  })

  output$dl <- downloadHandler(
    filename = function() paste0("columbus_", analysis()$method, "_", analysis()$var, ".csv"),
    content  = function(file) write.csv(analysis()$table, file, row.names = FALSE)
  )

  # ===========================================================================
  # REGRESSION
  # ===========================================================================
  reg <- eventReactive(input$reg_run, {
    validate(
      need(length(input$reg_x) >= 1, "Select at least one explanatory variable."),
      need(!(input$reg_y %in% input$reg_x),
           "The dependent variable cannot also be an explanatory variable.")
    )
    w <- build_lw(input$reg_nb, input$reg_k, "W")
    f <- reformulate(input$reg_x, response = input$reg_y)

    ols   <- lm(f, data = col_df)
    tests <- tryCatch(run_lm_tests(ols, w), error = function(e) NULL)
    slm   <- tryCatch(spatialreg::lagsarlm(f, data = col_df, listw = w, zero.policy = TRUE),
                      error = function(e) e)
    sem   <- tryCatch(spatialreg::errorsarlm(f, data = col_df, listw = w, zero.policy = TRUE),
                      error = function(e) e)

    list(formula = f, y = input$reg_y, x = input$reg_x, tests = tests,
         models = list(ols = summarise_model(ols, "ols", w),
                       slm = summarise_model(slm, "slm", w),
                       sem = summarise_model(sem, "sem", w)))
  }, ignoreNULL = FALSE)

  reg_sel <- reactive({
    r <- reg()
    m <- r$models[[input$reg_model]]
    validate(need(is.null(m$error),
                  paste("This model could not be estimated:", m$error)))
    m
  })

  output$reg_boxes <- renderUI({
    m <- reg_sel(); alpha <- input$reg_alpha
    ok <- m$mP >= alpha
    layout_columns(fill = FALSE, col_widths = c(3, 3, 3, 3),
      value_box("AIC", sprintf("%.2f", m$aic), showcase = icon("ruler"), theme = "primary"),
      value_box("Log-likelihood", sprintf("%.2f", m$ll), showcase = icon("chart-area"), theme = "secondary"),
      value_box("Residual Moran's I", sprintf("%.3f", m$mI),
                tags$small(paste("p =", fmt_p(m$mP))),
                showcase = icon(if (ok) "circle-check" else "triangle-exclamation"),
                theme = if (ok) "success" else "warning"),
      value_box(m$param_lab, sprintf("%.3f", m$param), showcase = icon("sliders"), theme = "info"))
  })

  output$reg_hyp_title <- renderText({
    paste("Hypotheses being tested \u2013", model_labels[[input$reg_model]])
  })

  output$reg_hyp <- renderUI({
    r <- reg(); m <- reg_sel(); alpha <- input$reg_alpha; k <- input$reg_model
    tab <- m$table
    is_sp <- grepl("^(rho|lambda)", tab$Term)
    sigx <- tab$Term[!is.na(tab$p) & tab$p < alpha & tab$Term != "(Intercept)" & !is_sp]
    coef_txt <- sprintf(
      "<b>Explanatory variables:</b> H<sub>0</sub>: \u03b2 = 0 for each of %s (no linear effect on %s) versus H<sub>1</sub>: \u03b2 \u2260 0. At \u03b1 = %.2f we reject H<sub>0</sub> for: <b>%s</b>.",
      paste(r$x, collapse = ", "), r$y, alpha,
      if (length(sigx)) paste(sigx, collapse = ", ") else "none of them")

    if (k == "ols") {
      sp_txt <- sprintf(
        "<b>Spatial independence of residuals (an OLS assumption):</b> H<sub>0</sub>: the residuals are spatially random (Moran's I = 0) versus H<sub>1</sub>: the residuals are spatially autocorrelated. Moran's I = %.3f, p = %s, so we %s.%s",
        m$mI, fmt_p(m$mP),
        if (m$mP < alpha) "<b>reject H<sub>0</sub></b>" else "<b>fail to reject H<sub>0</sub></b>",
        if (m$mP < alpha) " OLS assumptions are violated, so a spatial model is worth considering." else " OLS is not contradicted by the residual pattern.")
    } else {
      sp <- tab[is_sp, ][1, ]
      if (k == "slm") {
        sym <- "\u03c1"
        h0 <- sprintf("the value of %s in a neighbourhood is not affected by the values in neighbouring neighbourhoods (no spatial spillover)", r$y)
        h1 <- "neighbouring values do affect it (spatial lag dependence)"
      } else {
        sym <- "\u03bb"
        h0 <- "the error terms are not spatially correlated (no omitted, spatially structured factors)"
        h1 <- "the errors are spatially autocorrelated"
      }
      sp_txt <- sprintf(
        "<b>Spatial parameter:</b> H<sub>0</sub>: %s = 0, i.e. %s, versus H<sub>1</sub>: %s \u2260 0, i.e. %s. Estimate %s = %.3f (z = %.2f, p = %s), so we %s.",
        sym, h0, sym, h1, sym, sp$Estimate, sp$Stat, fmt_p(sp$p),
        if (!is.na(sp$p) && sp$p < alpha) "<b>reject H<sub>0</sub></b>" else "<b>fail to reject H<sub>0</sub></b>")
    }
    div(class = "interp", HTML(paste(coef_txt, sp_txt, sep = "<br><br>")))
  })

  output$reg_cmp_best <- renderUI({
    r <- reg(); alpha <- input$reg_alpha
    aics <- vapply(r$models, function(m) if (is.null(m$error)) m$aic else NA_real_, numeric(1))
    req(any(!is.na(aics)))
    best <- names(which.min(aics)); m <- r$models[[best]]
    srt <- sort(aics)
    gap <- if (length(srt) > 1) sprintf(", %.2f AIC points below the next-best model", srt[[2]] - srt[[1]]) else ""
    why <- switch(best,
      ols = "Adding a spatial term did not improve the fit enough, suggesting spatial dependence is weak once the explanatory variables are included.",
      slm = "The spatial lag term improves the fit, suggesting the outcome in a neighbourhood is influenced by the outcome in neighbouring areas (spillover).",
      sem = "The spatial error term improves the fit, suggesting spatially clustered factors not in the model also drive the outcome.")
    resid <- if (m$mP >= alpha) "Its residuals show no significant spatial autocorrelation."
             else "Its residuals still show some spatial autocorrelation."
    div(class = "best-model mt-2",
        HTML(sprintf("<b>\u2605 Best model: %s</b> (lowest AIC = %.2f%s). %s %s",
                     model_labels[[best]], m$aic, gap, why, resid)))
  })

  output$reg_coef_title <- renderText({
    paste("Coefficients \u2013", model_labels[[input$reg_model]])
  })

  output$reg_coef <- renderDT({
    m <- reg_sel()
    t <- m$table
    d <- data.frame(Term = t$Term,
                    Estimate = round(t$Estimate, 4),
                    `Std. error` = round(t$SE, 4),
                    `t / z value` = round(t$Stat, 3),
                    `p-value` = round(t$p, 4),
                    Sig = stars(t$p),
                    check.names = FALSE)
    datatable(d, rownames = FALSE, options = list(dom = "t", pageLength = 15, scrollX = TRUE))
  })

  output$reg_compare <- renderTable({
    r <- reg()
    rows <- lapply(names(r$models), function(k) {
      m <- r$models[[k]]
      if (!is.null(m$error)) {
        data.frame(Model = model_labels[[k]], AIC = NA_real_, dAIC = NA_real_,
                   LogLik = NA_real_, `Resid. Moran's I` = NA_real_, p = NA_character_,
                   check.names = FALSE)
      } else {
        data.frame(Model = model_labels[[k]], AIC = m$aic, dAIC = NA_real_,
                   LogLik = m$ll, `Resid. Moran's I` = m$mI, p = fmt_p(m$mP),
                   check.names = FALSE)
      }
    })
    d <- do.call(rbind, rows)
    d$dAIC <- d$AIC - min(d$AIC, na.rm = TRUE)
    d$AIC <- round(d$AIC, 2); d$dAIC <- round(d$dAIC, 2); d$LogLik <- round(d$LogLik, 2)
    d$`Resid. Moran's I` <- round(d$`Resid. Moran's I`, 3)
    if (any(!is.na(d$AIC))) {
      b <- which.min(d$AIC)
      d$Model[b] <- paste0(d$Model[b], " \u2605 best")
    }
    d
  }, striped = TRUE, width = "100%", na = "\u2013")

  output$reg_map <- renderPlot({
    m <- reg_sel()
    d <- columbus; d$res <- m$res
    lim <- max(abs(d$res))
    ggplot(d) +
      geom_sf(aes(fill = res), colour = "white", linewidth = 0.3) +
      scale_fill_gradient2(low = "#2c7bb6", mid = "#f7f7f7", high = "#d7191c",
                           midpoint = 0, limits = c(-lim, lim), name = "Residual") +
      guides(fill = guide_colourbar(barwidth = 14, barheight = 0.6, title.position = "top")) +
      labs(title = model_labels[[input$reg_model]],
           subtitle = "Red = under-predicted, blue = over-predicted") +
      coord_sf(datum = NA) +
      theme_map()
  })

  output$reg_scatter <- renderPlot({
    m <- reg_sel()
    w <- build_lw(input$reg_nb, input$reg_k, "W")
    sc <- quadrant_df(m$res, w)
    ggplot(sc, aes(z, lag)) +
      geom_hline(yintercept = 0, colour = "grey60") +
      geom_vline(xintercept = 0, colour = "grey60") +
      geom_point(colour = "#1d2b45", size = 3, alpha = 0.85) +
      geom_smooth(method = "lm", se = FALSE, colour = "#d7191c", linewidth = 1,
                  formula = y ~ x) +
      labs(x = "Residual (standardised)", y = "Spatial lag of residual (standardised)",
           subtitle = sprintf("Residual Moran's I = %.3f (p = %s)", m$mI, fmt_p(m$mP))) +
      theme_chart()
  })

  output$reg_lm <- renderTable({
    r <- reg()
    validate(need(!is.null(r$tests),
                  "Lagrange multiplier tests are not available in the installed spdep/spatialreg versions."))
    t <- r$tests
    data.frame(Test = t$Test, Statistic = round(t$Statistic, 3), df = t$df,
               `p-value` = fmt_p(t$p), Significant = ifelse(t$p < input$reg_alpha, "Yes", "No"),
               check.names = FALSE)
  }, striped = TRUE, width = "100%")

  output$reg_conclusion <- renderUI({
    r <- reg(); alpha <- input$reg_alpha
    aics <- vapply(r$models, function(m) if (is.null(m$error)) m$aic else NA_real_, numeric(1))
    req(any(!is.na(aics)))
    best <- names(which.min(aics)); m <- r$models[[best]]
    y_lab <- tolower(short_label(r$y))
    xs  <- m$table[m$table$Term %in% r$x, , drop = FALSE]
    sig <- !is.na(xs$p) & xs$p < alpha

    items <- vapply(seq_len(nrow(xs)), function(i) {
      lab <- short_label(xs$Term[i])
      if (sig[i]) {
        sprintf("<li><b>%s:</b> neighbourhoods with higher %s tend to have <b>%s %s</b>.</li>",
                lab, tolower(lab), if (xs$Estimate[i] < 0) "lower" else "higher", y_lab)
      } else {
        sprintf("<li><b>%s:</b> no clear link with %s once the other factor is taken into account.</li>",
                lab, y_lab)
      }
    }, character(1))

    spatial_txt <- switch(best,
      ols = sprintf("The simple model (OLS) turned out to be good enough. Once these factors are accounted for, there is no strong geographic pattern left over in %s, so a spatial model is not needed.", y_lab),
      slm = sprintf("The best model is the <b>spatial lag model</b>. In everyday terms: %s in one neighbourhood is linked to %s in the neighbourhoods around it, so it spills over across borders instead of staying inside one area.", y_lab, y_lab),
      sem = "The best model is the <b>spatial error model</b>. In everyday terms: something we did not measure (for example local policing, land use or the character of an area) affects neighbouring areas in a similar way, so the leftover differences are clustered rather than random.")
    resid_txt <- if (m$mP >= alpha) {
      "After using this model, the leftover differences between neighbourhoods look random, which is a good sign."
    } else {
      "Some geographic clustering is still left over, so the model does not capture everything."
    }

    short_txt <- if (any(sig)) {
      sprintf("%s differs across Columbus neighbourhoods in a way that is linked to %s%s.",
              cap(y_lab),
              paste(tolower(vapply(xs$Term[sig], short_label, character(1))), collapse = " and "),
              if (best != "ols") ", and it also clusters geographically, so neighbouring areas tend to look alike" else "")
    } else {
      sprintf("%s is not clearly explained by %s in these data%s.",
              cap(y_lab),
              paste(tolower(vapply(r$x, short_label, character(1))), collapse = " or "),
              if (best != "ols") ", although it does cluster geographically" else "")
    }

    div(class = "interp", HTML(paste0(
      "<p><b>What the explanatory factors tell us</b> (best model: ", model_labels[[best]], "):</p><ul>",
      paste(items, collapse = ""), "</ul>",
      "<p>", spatial_txt, " ", resid_txt, "</p>",
      "<p class='mb-0'><b>In short:</b> ", short_txt,
      " These results show association, not proof of cause and effect, and apply only to these 49 neighbourhoods.</p>")))
  })

  output$reg_text <- renderUI({
    r <- reg(); alpha <- input$reg_alpha
    parts <- character(0)
    ols <- r$models$ols
    if (is.null(ols$error)) {
      parts <- c(parts, sprintf(
        "Dependent variable: <b>%s</b>; explanatory variables: <b>%s</b>. OLS residuals show %s (Moran's I = %.3f, p = %s).",
        r$y, paste(r$x, collapse = ", "),
        if (ols$mP < alpha) "<b>significant spatial autocorrelation</b>, so OLS assumptions are violated"
        else "no significant spatial autocorrelation",
        ols$mI, fmt_p(ols$mP)))
    }
    if (!is.null(r$tests)) parts <- c(parts, decision_text(r$tests, alpha))
    aics <- vapply(r$models, function(m) if (is.null(m$error)) m$aic else NA_real_, numeric(1))
    if (any(!is.na(aics))) {
      best <- names(which.min(aics))
      parts <- c(parts, sprintf("By AIC, the best-fitting model is the <b>%s</b> (AIC = %.2f).",
                                model_labels[[best]], aics[[best]]))
    }
    sel <- r$models[[input$reg_model]]
    if (is.null(sel$error)) {
      sig_terms <- sel$table$Term[!is.na(sel$table$p) & sel$table$p < alpha &
                                    !(sel$table$Term %in% c("(Intercept)"))]
      parts <- c(parts, sprintf("In the displayed model (%s), significant terms at \u03b1 = %.2f: <b>%s</b>.",
                                model_labels[[input$reg_model]], alpha,
                                if (length(sig_terms)) paste(sig_terms, collapse = ", ") else "none"))
    }
    div(class = "interp", HTML(paste(parts, collapse = "<br><br>")))
  })

  # ===========================================================================
  # DATA
  # ===========================================================================
  output$data_tbl <- renderDT({
    datatable(col_df, rownames = FALSE,
              options = list(pageLength = 15, scrollX = TRUE))
  })
}

shinyApp(ui, server)

# Columbus
This dashboard works through three connected questions:

Global pattern: is crime spatially autocorrelated across the city? (Moran's I, Geary's C)
Local pattern: where are the clusters and spatial outliers? (LISA)
Explanation: does a spatial regression (lag or error) describe crime better than ordinary least squares?
Null hypothesis: values are distributed randomly in space (no spatial autocorrelation).

About the data
The classic Columbus, Ohio dataset (Anselin, 1988): 49 neighbourhoods (polygons) with crime and socio-economic indicators. Loaded from the spData package.

CRIME: residential burglaries and vehicle thefts per 1,000 households
HOVAL: housing value ($1,000s)
INC: household income ($1,000s)
OPEN: amount of open space in the neighbourhood
PLUMB: % of housing units without plumbing
DISCBD: distance to the central business district
Coordinates are planar and arbitrary (no real-world map projection), so maps show relative shape and position only.

Overview of the dashboard code

The whole dashboard is one file, app.R, in the standard Shiny pattern of a UI (what the user sees) and a server (the calculations behind it). The file runs from top to bottom in four blocks.

1. Data (back-end)

It loads the Columbus, Ohio shapefile (49 neighbourhoods) from the spData package as an sf spatial object.
The data has no real-world coordinate system, so the code strips the placeholder one. This stops ggplot2 trying to reproject it.
It also makes a rescaled copy of the polygons for the interactive Explore map only.
All statistics use the original geometry.

2. Helper functions

Build the neighbour structures (Queen, Rook or k-nearest) and the spatial weights.
Draw the neighbour-link map.
Compute Moran quadrants.
Format p-values.
Fit and summarise the regression models, including the Lagrange multiplier tests.

3. UI
It is a bslib page with four tabs, each with its own left-hand settings panel (frozen under a frozen top bar), and custom CSS for the navy-grey theme.

Explore: research question, data description, an interactive map with pop-ups, histogram, summary statistics, interpretations, and a downloadable data table.
Analysis: you pick a variable, neighbour definition and method, then press Calculate. You get a hypothesis card, value boxes, plots, a neighbourhood-structure map, an interpretation and a results table with CSV download.
Regression: you pick the variables and weights, then press Estimate models. It fits OLS, a spatial lag model (SLR) and a spatial error model (SEM). It shows hypotheses, coefficients, a model comparison with the best model marked, residual plots, diagnostics, an interpretation and a plain-English conclusion.
Data: the full attribute table.

4. Server

Settings are reactive inputs. The heavy calculations sit inside eventReactive() blocks triggered by the buttons, so nothing recomputes until you click.
The three methods are Global Moran’s I, Local Moran’s I (LISA) and Geary’s C.
p-values come from the normal approximation, and the weights are row-standardised.
The text boxes (hypotheses, interpretations, conclusion) are built with sprintf() from the live results, so they change with your choices.
Packages used
Package	What it does here
shiny	The reactive web-app framework
bslib	Page layout, cards, value boxes, theme and sidebars
sf	Spatial data handling, centroids, plotting polygons
spdep	Neighbours, spatial weights, Moran’s I, Geary’s C, LISA, LM tests
spData	Supplies the Columbus dataset
spatialreg	Estimates the spatial lag and spatial error models
ggplot2	Static plots (histogram, scatterplots, residual maps, z-plot)
leaflet	The interactive Explore map with hover labels and click pop-ups
DT	Searchable, sortable data tables
dplyr	Loaded for data handling, with light use
Deployment tools (not part of the app)
rsconnect creates manifest.json and publishes to shinyapps.io.
gert, with a GitHub token, pushes the files to GitHub.
These live in publish_all.R, not in the dashboard.

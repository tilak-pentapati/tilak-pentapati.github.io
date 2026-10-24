# Interactivity patterns for project pages

Working notes, not a page: `_quarto.yml` renders only `*.qmd`, so this file is
never published. Everything here works on static GitHub Pages (no server).

Reference implementations (Phase 4):

| Page | Patterns used |
|---|---|
| `r/StrokePredictionModel` | R model → `ojs_define()` → OJS calculator with in-browser verification |
| `python/satellites_overhead` | `panel-tabset` Map / Chart / Table, Python → `ojs_define()` → OJS chart, `itables`, light folium |
| `r/airbnb_sf` | `tmap_leaflet()` interactive maps with radio layers, Interactive / Static tabsets |
| `python/SP500_comp` | Plotly via `show_fig()`, `code-annotations`, `@fig-` / `@eq-` cross-references |

## 0. Ground rules

- **Freeze.** Editing a `.qmd` (front matter included) changes its freeze hash, so the
  page re-executes on the next render and needs its local data under `data/`. Render
  only the page you changed:
  - R: `QUARTO_R=/Users/jay/anaconda3/envs/r_env_gds/bin quarto render projects/r/<p>/index.qmd`
  - Python: `QUARTO_PYTHON=/Users/jay/anaconda3/envs/envs563/bin/python quarto render projects/python/<p>/index.qmd`
  - Commit the refreshed `_freeze/projects/.../` entry; the publish workflow has no R/Python.
- **Theme CSS.** Rules shared by these patterns live in `assets/styles/custom.scss`
  (DataTables container, inline code in annotation lists, KaTeX overflow, long URLs).
  Editing it changes the Bootstrap CSS hash, so do a full-site render afterwards (from
  freeze, nothing re-executes).
- **Weight.** Aim for < 1.5 MB HTML per page (GitHub Pages gzips it ~4x). Check with
  `wc -c docs/projects/<lang>/<p>/index.html` and `gzip -9 < … | wc -c`. Big offenders
  are embedded GeoJSON and inlined JS libraries.
- **Check.** Browser at ~375–500 px width, light and dark, `document.documentElement.scrollWidth <= innerWidth`,
  and no console errors. Serve `docs/` with
  `/Users/jay/anaconda3/envs/envs563/bin/python -m http.server 4321 --directory docs`.
- **Prefer** OJS (Observable Inputs + Plot, bundled with Quarto, no extra library) for
  calculators and filterable charts; Plotly when you want zoom/hover on a Python plot;
  Leaflet (tmap / folium) for maps. Shinylive only when the computation cannot be
  reduced to exported numbers (it downloads a WebAssembly runtime, tens of MB).

## 1. Export results to Observable JS (`ojs_define`)

R chunk (data frames arrive column-oriented; `transpose()` them in OJS):

```r
ojs_define(
  model_terms = model_terms,                    # data frame
  model_meta  = list(intercept = b0, threshold = p_star,
                     range = map(select(train, age, bmi), range))
)
```

Python chunk (`ojs_define` is injected by Quarto's Jupyter engine; no import):

```python
ojs_define(sat_long=sat_long, sat_times=times)  # pandas DataFrame + list
```

OJS: keep data/helper cells in their own chunk with `output: false`, otherwise every
`name = value` cell prints an inspector box:

````markdown
```{ojs}
//| code-summary: "Observable JS: model scoring"
//| output: false
terms = transpose(model_terms)
result = score(patient)
```

```{ojs}
//| code-summary: "Observable JS: inputs"
viewof patient = Inputs.form({
  age: Inputs.range([18, 82], {label: "Age (years)", step: 1, value: 65}),
  smoker: Inputs.radio(["No", "Yes"], {label: "Smoker", value: "No"}),
  work: Inputs.select(model_meta.levels.work_type, {label: "Work type"})
})
```
````

- `Inputs.*` render a real `<label for=…>`, so they are accessible as is. Give
  dynamic results `aria-live="polite"`.
- Bound sliders to the training range (export it) to avoid extrapolation.
- **Verify the port.** Export a sample of held-out rows with R/Python predictions and
  recompute them in OJS; print the max absolute difference on the page (the stroke page
  shows ≈1e-15). This is what makes the browser calculator trustworthy.
- For a fitted `tidymodels` workflow: `extract_recipe(fit)` then
  `tidy(prepped, number = i)` gives `step_dummy` levels and `step_normalize` means/SDs;
  `recipes::dummy_names(var, level)` reproduces the column names; join to
  `tidy(extract_fit_parsnip(fit))` and `stopifnot(!anyNA(...))` so a mismatch fails the render.
- Output cards: use Bootstrap classes that follow the theme (`card`, `alert-*`,
  `text-bg-*`, `text-body-secondary`), not hard-coded light colours.

OJS bar chart with filters (works in tabs, resizes with `width`):

```js
Plot.plot({
  width: Math.min(width, 820), height: 40 + 22 * rows.length, marginLeft: 170,
  x: {label: "Satellites overhead →", grid: true},
  y: {label: null, domain: rows.map(d => d.name),
      tickFormat: d => d.length > 24 ? d.slice(0, 23) + "…" : d},
  marks: [Plot.barX(rows, {x: "value", y: "name", fill: "#2e8b57", tip: true}), Plot.ruleX([0])]
})
```

Do not put `fig-cap` on an OJS chunk in pages whose ggplot titles hard-code
"Figure N:" — Quarto would number it and clash; write the caption as text instead.

## 2. Tabsets for alternatives

```markdown
::: {.panel-tabset}

## Map
(folium / leaflet chunk — keep the map in the first, visible tab so Leaflet sizes correctly)

## Chart
(OJS or Plotly)

## Table
(itables)

:::
```

Interactive / Static pairs use the same structure (`## Interactive`, `## Static`), so
the static figure stays available for print and for readers without JS.

## 3. Interactive tables (Python, `itables`)

`itables` is installed in `envs563` (2.x). Connected mode loads DataTables from a CDN
instead of inlining ~700 KB:

```python
from itables import show
show(df, connected=True, showIndex=False, scrollX=True,
     pageLength=10, classes="display compact", maxBytes=0)
```

`maxBytes=0` disables down-sampling (fine for a few hundred rows; keep tables small).
The dark-mode fix (transparent `.dt-container`) is already in `custom.scss`.
For R, the equivalent is `DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE))`
(check that `DT` is installed in `r_env_gds` first).

## 4. Plotly in Python

The notebook renderers add plotly.js **and MathJax** to every figure, and with Quarto's
require.js this produced console errors. Emit plain HTML instead, loading plotly.js once:

```python
import plotly.graph_objects as go
from plotly.offline import get_plotlyjs_version
from IPython.display import HTML, display

def style_fig(fig, height=420):
    """Transparent background, mid-grey text: readable in light and dark themes."""
    fig.update_layout(template="none", height=height, margin=dict(l=60, r=20, t=30, b=50),
                      paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                      font=dict(color="#8a8f98"), hovermode="x unified",
                      legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0))
    fig.update_xaxes(gridcolor="rgba(128,128,128,0.25)", zeroline=False)
    fig.update_yaxes(gridcolor="rgba(128,128,128,0.25)", zeroline=False)
    return fig

_plotly_loaded = False
def show_fig(fig):
    global _plotly_loaded
    html = fig.to_html(full_html=False, include_plotlyjs=False, include_mathjax=False,
                       config={"responsive": True, "displaylogo": False})
    if not _plotly_loaded:   # hide define() so plotly.js sets window.Plotly
        html = ('<script>window.__define = window.define; window.define = undefined;</script>'
                f'<script src="https://cdn.plot.ly/plotly-{get_plotlyjs_version()}.min.js"></script>'
                '<script>window.define = window.__define;</script>') + html
        _plotly_loaded = True
    display(HTML(html))
```

Pair with a cross-referenceable figure:

````markdown
```{python}
#| label: fig-scree
#| fig-cap: "Scree plot … Hover for exact values."
show_fig(style_fig(fig))
```
````

Swapping the plotted series without re-rendering: `updatemenus` with
`method="update"` (see `fig-loadings` on the SP500 page).

## 5. Interactive maps in R (tmap 4 view mode)

```r
# helpers (once per page)
web_layer <- function(x, ...) {
  x <- x %>% select(...) %>%
    st_simplify(preserveTopology = TRUE, dTolerance = 200) %>%   # units of the projected CRS
    st_transform(4326)
  g <- st_geometry(x)                                            # round to ~1 m
  st_geometry(x) <- st_sfc(lapply(g, function(p)
    structure(rapply(unclass(p), round, how = "replace", digits = 5), class = class(p))),
    crs = st_crs(g))
  x
}
web_basemap <- tm_basemap("Esri.WorldGrayCanvas", group.control = "none")  # keyless; stays on
web_view <- tm_view(set_view = c(lon, lat, zoom), control.collapse = FALSE)
```

````markdown
```{r map_web}
#| out-width: "100%"
#| out-height: "480px"
lyr <- x %>% mutate(label = sprintf("ZIP %s: %s", ZIP, formatC(value, format = "d", big.mark = ","))) %>%
  web_layer(label, value, other)
m <- web_basemap + tm_shape(lyr) +
  tm_polygons(fill = "value", fill.scale = tm_scale_intervals(style = "fisher", n = 5, values = "brewer.yl_gn"),
              fill.legend = tm_legend("Value"), fill_alpha = 0.8, col = "white", lwd = 0.3,
              id = "label", popup = FALSE, group = "Value", group.control = "radio") +
  tm_polygons(fill = "other", ..., group = "Other", group.control = "radio") +
  web_view
m_lf <- tmap_leaflet(m, show = FALSE)   # returns invisibly: assign, then print
m_lf
```
````

Gotchas found on `airbnb_sf`:

- `tmap_leaflet()` returns invisibly, and `x %>% identity()` was still dropped by knitr
  on that page; assigning and printing the object works.
- Widget size comes from `out-width` / `out-height` chunk options; setting `$height` on
  the widget is ignored.
- With `group.control = "radio"`, the basemap must have `group.control = "none"`, or
  choosing a layer switches the basemap off. Legends follow their layer automatically.
- Every layer embeds its own copy of the geometry, and tmap popups add about 700 bytes
  per polygon. Hover labels (`id = "label"`) plus simplification and rounding took the
  page from 3.0 MB to 1.3 MB.
- Do not use `format()` inside `sprintf()` for labels (it pads to a common width); use
  `formatC(x, format = "d", big.mark = ",")`.
- Keep the static `tmap_mode("plot")` maps for the Static tab; `tmap_leaflet()` does not
  change the global mode.

## 6. Interactive maps in Python (folium)

```python
map_gdf = gdf[['name', 'value', 'geometry']].copy()          # displayed columns only
map_gdf['geometry'] = shapely.set_precision(
    map_gdf.geometry.simplify(0.05, preserve_topology=True).values, 0.001)

m = folium.Map(location=[35, 0], zoom_start=2,
               tiles='https://server.arcgisonline.com/ArcGIS/rest/services/Canvas/World_Light_Gray_Base/MapServer/tile/{z}/{y}/{x}',
               attr='Tiles &copy; Esri')
ch = folium.Choropleth(geo_data=map_gdf, data=map_gdf, columns=('name', 'value'),
                       key_on='feature.properties.name', fill_color='Greens', highlight=True)
ch.geojson.add_child(folium.features.GeoJsonTooltip(fields=['name']))   # same layer, no 2nd copy
ch.geojson.add_child(folium.features.GeoJsonPopup(fields=['name', 'value']))
ch.add_to(m)
m
```

Folium maps sit in an `<iframe srcdoc>`, so the GeoJSON is HTML-escaped (about 2x).
On `satellites_overhead` a single layer, 0.05° simplification and a 0.001° grid took the
page from 18.4 MB to 1.3 MB. CARTO tiles now need an API key; Mapbox needs a token, so
use the keyless Esri Canvas tiles (as on the home page).

## 7. Cross-references and code annotations

- Figures: `#| label: fig-xyz` + `#| fig-cap:` on the chunk, then `@fig-xyz` in the text.
  Labels must start with `fig-`; this works for Plotly HTML output too.
- Equations: `$$ … $$ {#eq-xyz}`, then `@eq-xyz`. Display maths inside a list item also
  works (indent it under the item). `.katex-display` scrolls sideways on phones (custom.scss).
- Tables: `#| label: tbl-xyz` + `#| tbl-cap:` (R `kable`, Python DataFrame output).
- Code annotations: `code-annotations: below` in the page front matter, `# <1>` at the
  end of a code line, then a numbered list straight after the chunk:

````markdown
```{python}
pca = PCA()              # <1>
pca.fit(returns_table)   # <2>
```
1. No `n_components`: all components are kept.
2. Mean-centres and takes the SVD; see @eq-eigen.
````

Wide markdown tables: wrap them in `::: {.table-responsive}` … `:::`.

## 8. Suggestions for the remaining projects

| Project | Suggested upgrade |
|---|---|
| `r/GCH_census` (OLS, Green Central Heating) | OJS "what-if" predictor from `coef()` + model-matrix means (pattern 1); `tbl-` cross-refs for the 10 `kable` tables; wrap wide tables in `.table-responsive` (known 375 px overflow). |
| `r/aps_eda` (SOC × address duration) | Plotly/ggplotly or OJS filter on the main bivariate plots; `.table-responsive` for the wide tables (known overflow). |
| `r/se_m_airbnb` (Madrid point patterns, spatial econometrics) | `tmap_leaflet()` map of listings/intensity with radio layers (pattern 5) in an Interactive / Static tabset; `@eq-` labels on the SAR/SEM equations. |
| `r/stm_covid` (multilevel, GWR, spatio-temporal) | Interactive GWR coefficient map, one radio layer per coefficient (pattern 5); OJS time slider over the weekly rates (pattern 1). |
| `python/bd_gps` (GeoLife, PySpark) | Pre-aggregate trajectories to a small grid, then a folium/OJS heat map; never ship raw points (weight). |
| `python/retail_site_select` (catchments) | Folium map with catchment and candidate-site layers (pattern 6); OJS sliders for the scoring weights so readers can re-rank sites. |
| `python/stPopLCR` (population trends) | Map / Chart / Table tabset: choropleth by year, OJS line chart with a ward filter, `itables` of the trend table (patterns 2–3). |

![](readme/fig1.png)

[![CRAN status](https://www.r-pkg.org/badges/version/rGEDI)](https://CRAN.R-project.org/package=rGEDI)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/rGEDI)](https://cran.r-project.org/package=rGEDI)
![Licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg)

# rGEDI: NASA's Global Ecosystem Dynamics Investigation (GEDI) Data Visualization and Processing

**Authors:** Carlos Alberto Silva and Caio Hamamura.

`rGEDI` searches, downloads, streams, reads, filters, clips, summarizes,
simulates, models and maps GEDI data. The package supports Level 1B waveforms,
Level 2A/2B vegetation structure, Level 3 grids and Level 4A/4B biomass,
including integration with other sensors for wall-to-wall mapping.

## From laser measurements to maps

GEDI is a laser instrument aboard the International Space Station. Each
measurement samples a footprint approximately 25 m across. Its waveform
records returned laser energy through the vegetation and down to the ground.
GEDI samples tracks; it does not measure every location in a landscape.

This tutorial follows those measurements from **Level 1B waveforms**, through
**Level 2A/2B vegetation metrics**, **Level 3 grids**, and **Level 4A/4B biomass**.
We then simulate GEDI waveforms from measured ALS data and combine GEDI with
other satellite data to make a wall-to-wall map.

## Installation and configuration

### Install rGEDI

The released version is available from CRAN:

```r
install.packages("rGEDI")
```

The examples below use the **0.6.1 development API**. Install the GitHub version
to follow them:

```r
install.packages("remotes")
remotes::install_github("carlos-alberto-silva/rGEDI")

# Alternative: the development build from r-universe, once available.
# install.packages("rGEDI", repos = c(
#   "https://carlos-alberto-silva.r-universe.dev",
#   "https://cloud.r-project.org"))

library(rGEDI)
library(data.table)
library(sf)
library(terra)

outdir <- file.path(tempdir(), "rGEDI-tutorial")
dir.create(outdir, showWarnings = FALSE)
```

Bundled data support offline extraction, statistics and static plots.
Install the optional R packages needed for interactive imagery, measured ALS,
model selection and the figure layouts:

```r
install.packages(c("leaflet", "lidR", "plot3D", "randomForest",
                   "patchwork", "viridisLite", "rnaturalearthdata"))
```

Interactive imagery tiles require internet access. Local HDF5 and raster
workflows do not require Python or cloud authentication.

### Configure optional Python backends

The cloud tools use `reticulate` with `earthaccess` for Earthdata access,
`h5py` for streamed HDF5, and `earthengine-api` for Google Earth Engine.

```r
install.packages("reticulate")

# Optional: select an existing Python executable before loading reticulate.
# Sys.setenv(RETICULATE_PYTHON = "/path/to/your/python")
rGEDI_configure(install = TRUE)

# Include the optional gediDB backend when needed:
# rGEDI_configure(install = TRUE, gedidb = TRUE)
```

Inspect the selected interpreter and module availability:

```r
reticulate::py_config()
rGEDI_configure(install = FALSE)
gee_available()
# gedidb_available()  # If the optional database backend was configured.
```

### Authenticate Earthdata and Earth Engine

Earthdata access can use a standard `~/.netrc`, the `NETRC` environment
variable, or an interactive login. Keep credentials outside scripts.

```r
# Optional alternative to the standard ~/.netrc location:
# Sys.setenv(NETRC = "/secure/path/to/.netrc")
earthdata_login()
```

For GEE, supply a Google Cloud project registered for Earth Engine and follow
the authentication prompt on first use:

```r
Sys.setenv(EE_PROJECT = "your-google-cloud-project")
ee <- ee_initialize(project = Sys.getenv("EE_PROJECT"))
```

These setup calls are explicit, optional steps for cloud examples; package
installation and checks do not install Python modules or authenticate services.

## 1 Start with Level 1B waveforms

Level 1B contains geolocated waveforms. First open the small Cerrado example
distributed with the package.

```r
level1b_path <- unzip(system.file("extdata",
  "GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.zip",
  package = "rGEDI"), exdir = outdir)
level1b <- readLevel1B(level1b_path)

level1b_geo <- getLevel1BGeo(level1b)
head(level1b_geo[, .(shot_number, longitude_bin0, latitude_bin0)], 3)
```

```text
         shot_number longitude_bin0 latitude_bin0
               <i64>          <num>         <num>
1: 19640119100108615      -44.13999     -13.72638
2: 19640119300108616      -44.13967     -13.72596
3: 19640119500108617      -44.13936     -13.72555
```

### Plot the footprints

Coordinates tell us where the laser sampled. Plot them over imagery:

```r
footprint_map <- leaflet::leaflet(data = as.data.frame(level1b_geo)) |>
  leaflet::addProviderTiles("Esri.WorldImagery") |>
  leaflet::addCircleMarkers(
    lng = ~longitude_bin0, lat = ~latitude_bin0,
    radius = 2, color = "red", stroke = FALSE, fillOpacity = 1
  ) |>
  leaflet::addScaleBar(position = "topright",
    options = leaflet::scaleBarOptions(imperial = FALSE)) |>
  leaflet::addLegend(position = "topright", colors = "red",
    labels = "Samples", title = "GEDI Level 1B")
footprint_map
```

<p align="center"><img src="readme/fig2.PNG" width="800" alt="GEDI Level 1B footprints over imagery"></p>

### Examine one waveform

The shot identifier connects a waveform to its footprint and higher-level
metrics. Keep shot identifiers as character or integer64 values.

```r
shot <- level1b_geo$shot_number[1]
waveform <- getLevel1BWF(level1b, shot_number = shot)
plot(waveform)
```

<p align="center"><img src="readme/fig3.png" width="750" alt="GEDI waveform in raw and relative amplitude units"></p>

### Find additional granules

For your own area, search with a bounding box and date interval. Earthdata
downloads require an Earthdata account.

```r
study_area <- st_read(system.file("extdata", "stands_cerrado.shp",
                                  package = "rGEDI"), quiet = TRUE)
box <- st_bbox(study_area)
xmin <- unname(box["xmin"]); xmax <- unname(box["xmax"])
ymin <- unname(box["ymin"]); ymax <- unname(box["ymax"])

urls <- gedifinder("GEDI01_B", ul_lat = ymax, ul_lon = xmin,
  lr_lat = ymin, lr_lon = xmax, version = "002",
  daterange = c("2019-04-18", "2019-04-19"))
# earthdata_login()
# gediDownload(urls, outdir = outdir)
```

`gedifinder(orbit = "O01964", product = "GEDI02_A", version = "002")`
finds the four production granules for that orbit. `getTrack()` and
`plot_gedi_orbit_animation()` can display their time-ordered tracks.

<p align="center"><img src="readme/gedi-orbit-animation.gif" width="650" alt="GEDI sampling tracks from the International Space Station"></p>

## 2 Describe vegetation with Level 2A and Level 2B

Level 2A supplies ground elevation and relative-height metrics. RH98 is the
height associated with 98% of cumulative waveform energy, relative to the
estimated ground. Level 2B describes canopy cover, plant area index (PAI), and
plant area volume density (PAVD).

```r
level2a_path <- unzip(system.file("extdata",
  "GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
  package = "rGEDI"), exdir = outdir)
level2b_path <- unzip(system.file("extdata",
  "GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.zip",
  package = "rGEDI"), exdir = outdir)
level2a <- readLevel2A(level2a_path)
level2b <- readLevel2B(level2b_path)
level2a_metrics <- getLevel2AM(level2a)
level2b_vpm <- getLevel2BVPM(level2b)
head(level2a_metrics[, .(beam, shot_number, rh50, rh98, rh100)], 3)
```

```text
       beam       shot_number  rh50  rh98 rh100
     <char>             <i64> <num> <num> <num>
1: BEAM0001 19640119100108615 -0.14  3.25  4.30
2: BEAM0001 19640119300108616 -0.26  3.07  4.08
3: BEAM0001 19640119500108617 -0.26  3.22  4.41
```

```r
head(level2b_vpm[, .(beam, shot_number, pai, cover, fhd_normal)], 3)
```

```text
    beam       shot_number         pai       cover fhd_normal
   <int>             <i64>       <num>       <num>      <num>
1:     1 19640119100108615 0.056964017 0.028080121  1.3729203
2:     1 19640119300108616 0.003395777 0.001696448  0.6365142
3:     1 19640119500108617 0.037773777 0.018709604  1.4299893
```

### Connect height metrics to the waveform

```r
plotWFMetrics(level1b, level2a, shot_number = shot,
              rh = c(25, 50, 75, 90, 98))
```

<p align="center"><img src="readme/fig8.png" width="650" alt="Waveform, cumulative energy and relative-height metrics"></p>

### Explore vertical vegetation profiles

PAI describes plant area per ground area; PAVD describes its vertical
distribution. Extract the profiles before plotting a selected beam.

```r
pai_profile <- getLevel2BPAIProfile(level2b)
pavd_profile <- getLevel2BPAVDProfile(level2b)
profile_beam <- unique(pai_profile$beam)[1]
plotPAIProfile(pai_profile, beam = profile_beam)
plotPAVDProfile(pavd_profile, beam = profile_beam)
```

<p align="center"><img src="readme/fig9.png" width="1000" alt="PAI and PAVD profiles along GEDI tracks"></p>

### Filter and clip the footprints

Choose quality criteria appropriate to the analysis. This local example
retains high-quality, undegraded nighttime measurements from power beams.

```r
power_beams <- c("BEAM0101", "BEAM0110", "BEAM1000", "BEAM1011")
level2a_good <- level2a_metrics[
  quality_flag == 1 & degrade_flag == 0 &
  sensitivity >= 0.9 & sensitivity <= 1 &
  solar_elevation < 0 & beam %in% power_beams
]
level2b_good <- level2b_vpm[l2b_quality_flag == 1 & sensitivity >= 0.9]

# Extracted tables return clipped tables.
bbox <- c(xmin, xmax, ymin, ymax)
level2a_bbox <- clip(level2a_metrics, bbox)
level2a_inside <- clip(level2a_metrics, study_area)

# Open HDF5 objects can instead be written as subset HDF5 files.
level2a_h5_clip <- clip(level2a, study_area,
  output = file.path(outdir, "level2a-clipped.h5"))
```

<p align="center"><img src="readme/fig4.png" width="1000" alt="GEDI clipping by bounding box and stand polygons"></p>

### Summarize each polygon

Use a polygon identifier to keep the summaries separate. Here cover is a
fraction, not a percentage.

```r
polygons <- study_area
polygons$unit_id <- seq_len(nrow(polygons))
footprints <- clip(level2b_good, polygons, split_by = "unit_id")
structure_metrics <- function(pai, cover) list(
  n = length(pai), mean_pai = mean(pai), mean_cover = mean(cover)
)
polygon_summary <- polyStatsLevel2BVPM(
  footprints, ~structure_metrics(pai, cover), id = "poly_id"
)
polygon_summary
```

```text
   poly_id     n   mean_pai mean_cover
     <int> <int>      <num>      <num>
1:       2    98 0.11910636 0.05553518
2:       1   100 0.15160601 0.06934074
3:       5    52 0.06078231 0.02914990
4:       4    20 0.25203198 0.11343149
5:       3    29 0.23915693 0.10965394
```

### Build your own footprint grids

`gridStatsLevel2AM()` and `gridStatsLevel2BVPM()` calculate user-defined
statistics. Their geographic resolution is in degrees.

```r
grid_summary <- function(x) list(
  min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE),
  mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE)
)
rh100_grid <- gridStatsLevel2AM(level2a_metrics,
  func = grid_summary(rh100), res = 0.005)
pai_grid <- gridStatsLevel2BVPM(level2b_vpm,
  func = grid_summary(pai), res = 0.005)
```

<p align="center"><img src="readme/fig5.png" width="390" alt="RH100 grid statistics"><img src="readme/fig6.png" width="390" alt="PAI grid statistics"></p>

For an exact 1 km cell size, supply a projected template to `rasterize()`.
Choose a projection appropriate to your area; UTM 23S fits this example.

```r
points_utm <- terra::project(to_vect(level2a_metrics), "EPSG:32723")
bounds <- as.vector(terra::ext(points_utm))
bounds[c(1, 3)] <- floor(bounds[c(1, 3)] / 1000) * 1000
bounds[c(2, 4)] <- ceiling(bounds[c(2, 4)] / 1000) * 1000
template <- terra::rast(terra::ext(bounds), resolution = 1000,
                        crs = "EPSG:32723")
rh98_1km <- rGEDI::rasterize(level2a_metrics, "rh98",
  template = template, fun = "mean")
pai_1km <- rGEDI::rasterize(level2b_vpm, "pai",
  template = template, fun = "mean")
```

<p align="center"><img src="readme/fig-package-spatial.png" width="1100" alt="Polygon clipping and projected 1 km mean RH98 and PAI grids"></p>

These grids summarize sampled footprints. White cells have no included
measurements; they are not zero-height vegetation.

## 3 Read and explore native Level 3 grids

NASA also distributes nominal 1 km Level 3 summaries. These are different
from the user-defined grids above. A small regional subset for Sequoia National Forest,
California, is bundled with its USFS administrative boundary so this example runs offline.

```r
level3_file <- system.file("extdata", "GEDI03_Sequoia.tif", package = "rGEDI")
level4b_file <- system.file("extdata", "GEDI04_B_Sequoia.tif", package = "rGEDI")
level3 <- readLevel3(level3_file)
level4b <- readLevel4B(level4b_file)
```

```text
    product  rows columns layers resolution_m
     <char> <num>   <num>  <num>        <num>
1:  Level 3   178     146      2     1000.895
2: Level 4B   178     146      2     1000.895
```

### Clip a grid to a mapped boundary

Transform the boundary to the raster's CRS before using its coordinates.
`clipLevel3()` accepts either an extent or a polygon.

```r
grid_aoi <- st_read(system.file("extdata", "sequoia.shp",
                                package = "rGEDI"), quiet = TRUE)
native_aoi <- st_transform(grid_aoi, terra::crs(level3))
grid_bbox <- as.vector(terra::ext(terra::vect(native_aoi)))
level3_bbox <- clipLevel3(level3, grid_bbox)
level3_aoi <- clipLevel3(level3, terra::vect(native_aoi))

plotLevel3(level3[["rh100_mean"]], main = "Level 3 mean RH100 (m)")
plotLevel3(level3_aoi[["rh100_mean"]], main = "Mean RH100 inside Sequoia National Forest")
```

<details>
<summary>Complete map layout code: white and grey background, shared palette, boundary, north arrows and scale bars</summary>

```r
library(ggplot2)
library(patchwork)
outline <- st_transform(grid_aoi,32611)
draw <- function(full,clipped,product,unit,palette) {
  projector <- if(product=="3") projectLevel3 else projectLevel4B
  full <- projector(full,"EPSG:32611",method="near")
  clipped <- projector(clipped,"EPSG:32611",method="near")
  e <- as.vector(ext(full));dx<-e[2]-e[1];dy<-e[4]-e[3]
  limits <- range(values(full),na.rm=TRUE)
  panel <- function(r,title) {
    d <- as.data.frame(r,xy=TRUE,na.rm=TRUE);names(d)[3]<-"value"
    ggplot(d,aes(x,y,fill=value))+geom_tile(width=res(r)[1],height=res(r)[2])+
      geom_sf(data=outline,inherit.aes=FALSE,fill=NA,colour="white",linewidth=1.1)+
      geom_sf(data=outline,inherit.aes=FALSE,fill=NA,colour="#50575A",linewidth=.45)+
      scale_fill_viridis_c(option=palette,limits=limits,name=unit,
        guide=guide_colourbar(barwidth=grid::unit(7,"cm"),barheight=grid::unit(.35,"cm")))+
      coord_sf(crs=32611,xlim=e[1:2],ylim=c(e[3]-.12*dy,e[4]),datum=NA,expand=FALSE)+
      annotate("segment",x=e[1]+.08*dx,xend=e[1]+.08*dx+20000,y=e[3]-.07*dy,yend=e[3]-.07*dy,linewidth=.8)+
      annotate("text",x=e[1]+.08*dx+10000,y=e[3]-.035*dy,label="20 km",size=3.5)+
      annotate("rect",xmin=e[2]-.13*dx,xmax=e[2]-.03*dx,ymin=e[4]-.20*dy,ymax=e[4]-.015*dy,fill="white",colour="#E1E4E0",linewidth=.2)+
      annotate("segment",x=e[2]-.08*dx,xend=e[2]-.08*dx,y=e[4]-.18*dy,yend=e[4]-.07*dy,
        arrow=grid::arrow(type="closed",length=grid::unit(2,"mm")))+
      annotate("text",x=e[2]-.08*dx,y=e[4]-.035*dy,label="N",size=4)+
      labs(title=title,subtitle="Sequoia National Forest")+
      theme_void(base_size=13)+theme(plot.title=element_text(face="bold"),legend.position="bottom",
        panel.background=element_rect(fill="#F2F3F1",colour=NA),plot.background=element_rect(fill="white",colour=NA))
  }
  panel(full,paste0("A  Native Level ",product))+panel(clipped,"B  Clipped to the forest")+
    plot_layout(guides="collect")+plot_annotation(caption="Grey outline: USFS administrative boundary. Light grey: no displayed estimate. Native nominal 1 km cells; nearest-neighbour projection for display.") &
    theme(legend.position="bottom",plot.caption=element_text(size=9,hjust=0))
}

structure_map <- draw(level3[["rh100_mean"]], level3_aoi[["rh100_mean"]],
                      "3", "RH100 (m)", "viridis")
structure_map
ggsave(file.path(outdir,"level3-sequoia.png"),structure_map,
       width=11,height=7,dpi=250,bg="white")
```

</details>

<p align="center"><img src="readme/fig-level3-products.png" width="1050" alt="Native Level 3 mean RH100 before and after clipping to the Sequoia National Forest boundary"></p>

### Extract pixel values

`getLevel3()` returns the native cell measurements intersecting the polygon.
`ID` identifies the input polygon; each row is an intersected cell.

```r
level3_values <- getLevel3(level3, geometry = terra::vect(native_aoi))
head(level3_values[is.finite(rh100_mean)], 3)
```

```text
      ID rh100_mean counts
   <num>      <num>  <num>
1:     1   13.43818     11
2:     1   11.44869     61
3:     1   19.22667     18
```

## 4 Move from structure to biomass with Level 4A and Level 4B

Level 4A estimates aboveground biomass density (AGBD) at footprints.
Level 4B provides gridded estimates and uncertainty information. Both use
Mg/ha, equivalent to metric tonnes per hectare.

### Level 4A footprints through Earth Engine

This optional online step introduces the new mapping area, Shenandoah
National Park. Its NPS boundary is supplied as a shapefile.

```r
# Configure Python dependencies once, if needed:
# rGEDI_configure(install = TRUE)
# Sys.setenv(EE_PROJECT = "your-authenticated-project")
ee <- ee_initialize(project = Sys.getenv("EE_PROJECT"))
aoi <- st_read(system.file("extdata", "shenandoah.shp", package = "rGEDI"),
               quiet = TRUE)
box <- st_bbox(aoi)
address <- gedifinder("GEDI04_A", ul_lat = box[["ymax"]],
  ul_lon = box[["xmin"]], lr_lat = box[["ymin"]], lr_lon = box[["xmax"]],
  version = "002.1", daterange = c("2019-05-01", "2019-10-01"), source = "gee")
level4a <- readLevel4A(address)
lazy <- getLevel4A(level4a, quality = TRUE, collect = FALSE)
lazy <- clip(lazy, aoi)  # Apply the polygon filter before transferring data.
biomass_footprints <- getLevel4A(lazy, quality = TRUE, max_features = 10000,
  cols = c("beam", "shot_number", "agbd", "agbd_se", "sensitivity",
           "solar_elevation", "l4_quality_flag", "degrade_flag"))
head(biomass_footprints[, .(beam, shot_number, agbd, agbd_se, sensitivity)], 3)
```

```text
    beam       shot_number     agbd  agbd_se sensitivity
   <int>            <char>    <num>    <num>       <num>
1:    11 43081100300231185 453.4523 10.32448    0.977248
2:    11 43081100300231211 366.5885 10.31399    0.979268
3:    11 43081100300231207 398.6888 10.31586    0.989330
```

The finder-reader-extractor sequence also supports Earthdata Cloud and
gediDB. Select `source = "earthdata"` or `source = "gedidb"` and the product
version provided by that service. GEE retains numeric beam identifiers;
its power beams are `5`, `6`, `8`, and `11`.

### Plot and clip Level 4B grids offline

Return to the bundled Sequoia grids to compare regional biomass estimates
with the Level 3 structure metrics.

```r
level4b_bbox <- clipLevel4B(level4b, grid_bbox)
level4b_aoi <- clipLevel4B(level4b, terra::vect(native_aoi))
plotLevel4B(level4b[["MU"]], main = "Level 4B AGBD (Mg/ha)")
plotLevel4B(level4b_aoi[["MU"]], main = "AGBD inside Sequoia National Forest")
```

The same layout uses magma for biomass:

```r
biomass_map <- draw(level4b[["MU"]], level4b_aoi[["MU"]],
                    "4B", "AGBD (Mg/ha)", "magma")
biomass_map
ggsave(file.path(outdir,"level4b-sequoia.png"),biomass_map,
       width=11,height=7,dpi=250,bg="white")
```

<p align="center"><img src="readme/fig-level4b-products.png" width="1050" alt="Native Level 4B AGBD before and after polygon clipping"></p>

### Extract cells and calculate a polygon summary

```r
level4b_values <- getLevel4B(level4b, geometry = terra::vect(native_aoi))
head(level4b_values[is.finite(MU)], 3)
```

```text
      ID       MU        SE
   <num>    <num>     <num>
1:     1 57.40922  8.426760
2:     1 25.16608  5.206754
3:     1 54.39237 22.470741
```

```r
level4b_stats <- polyStats(level4b, terra::vect(native_aoi),
                           fun = mean, na.rm = TRUE)
level4b_stats
```

```text
      ID       MU       SE
   <int>    <num>    <num>
1:     1 92.57902 12.80514
```

`MU` is the mean AGBD estimate and `SE` is its cell-level standard error.
Here the summary averages the intersected cell values. The mean of the SE
values is not the uncertainty of the regional mean.

The corresponding clipping calls return:

```text
        output cells valid
        <char> <num> <num>
1:     L3 bbox 19908 19515
2:  L3 polygon 19908  4952
3:    L4B bbox 19908 19430
4: L4B polygon 19908  4942
```

## 5 Simulate GEDI waveforms from measured ALS

ALS resolves individual laser returns. `gediWFSimulator()` combines those
returns into a GEDI-like waveform at a chosen footprint location. The input
below is the original measured Amazon ALS example, not a random point cloud.
The package also includes `Savanna.laz`.

```r
als_file <- system.file("extdata", "Amazon.laz", package = "rGEDI")
als <- lidR::readLAS(als_file)
centre <- c(mean(range(als@data$X)), mean(range(als@data$Y)))

sim_clean <- gediWFSimulator(als, output = file.path(outdir, "sim-clean.h5"),
                            coords = centre, noise = 0, seed = 11)
sim_noisy <- gediWFSimulator(als, output = file.path(outdir, "sim-noisy.h5"),
                            coords = centre, noise = 0.001, seed = 11)
metrics_clean <- gediWFMetrics(sim_clean, rhRes = 1, ground_method = "lowest_peak")
metrics_noisy <- gediWFMetrics(sim_noisy, rhRes = 1, ground_method = "lowest_peak")
sim_waveform <- getLevel1BWF(sim_clean, metrics_clean$shot_number[1])
plot(sim_waveform, relative = TRUE)
```

<p align="center"><img src="readme/fig-package-als.png" width="1000" alt="Measured Amazon and Savanna ALS clouds and GEDI-like simulated waveforms"></p>

Repeating these calls for both measured clouds gives:

```text
      site noise     cover       rh50      rh90      rh98 rh100 waveEnergy
    <char> <num>     <num>      <num>     <num>     <num> <num>      <num>
1:  Amazon 0.000 0.7795247 20.0403346 30.525283 31.963053 36.15   1.000000
2:  Amazon 0.001 0.7767196 20.0361180 30.535343 31.995403 70.95   1.002188
3: Savanna 0.000 0.1565774  0.4287428  4.700784  7.047970 13.05   1.000000
4: Savanna 0.001 0.1549919  0.4276787  4.759650  7.226346 47.85   1.007346
```

Waveform heights use the estimated ground; ALS display heights use each
cloud's minimum elevation. Notice how noise affects RH100 more than RH98.
Ground estimation and the waveform noise window matter when interpreting
these metrics.

## 6 Integrate GEDI and other sensors for wall-to-wall mapping

GEDI footprints provide training observations, while spatially continuous
satellite and terrain predictors allow a model to estimate unsampled cells.
This example maps biomass in Shenandoah National Park using annual AlphaEarth
satellite embeddings and NASADEM terrain variables accessed through rGEDI.

### See global GEDI data through GEE

A small display preview of the native Level 4B product illustrates global
coverage without transferring the footprint archive.

```r
global_address <- gedifinder("GEDI04_B", ul_lat = 52, ul_lon = -180,
  lr_lat = -52, lr_lon = 180, version = "002", source = "gee")
global_biomass <- readLevel4B(global_address, metrics = "MU")
map_view(list(AGBD = global_biomass),
  vis = list(AGBD = list(min = 0, max = 400,
    palette = c("#000004", "#721F81", "#F1605D", "#FCFDBF"))))
```

<p align="center"><img src="readme/fig-gee-global.png" width="1100" alt="Global native GEDI Level 4B biomass viewed through Earth Engine"></p>

The figure uses a 20 km nearest-neighbour display preview. It is not a newly
estimated global biomass raster. The underlying product has its own sampling
support, period and uncertainty.

### Complete multi-year upscaling workflow

The code below runs directly, without sourcing another script. It retrieves
May–September GEDI observations from 2019, 2021 and 2022, extracts year-matched
AlphaEarth and NASADEM predictors, and holds out 2 km blocks across all years.
Use a new cache folder when changing the area, dates or selection settings.
The cache path uses the same code on Windows, macOS and Linux.

### Prepare the cache and global display

```r
root <- Sys.getenv("RGEDI_EXAMPLE_CACHE", file.path(path.expand("~"), "rGEDI-cache"))
dir.create(root, recursive = TRUE, showWarnings = FALSE)
dir.create("readme/results", recursive = TRUE, showWarnings = FALSE)
log <- function(...) {cat(format(Sys.time(), "%H:%M:%S"), ..., "\n"); flush.console()}
if (!nzchar(Sys.getenv("EE_PROJECT"))) stop("Set EE_PROJECT to an authenticated Earth Engine project.")
ee <- ee_initialize(project = Sys.getenv("EE_PROJECT"), authenticate = FALSE)
ee$data$setDeadline(180000L)
log("Earth Engine initialized")

# Native global Level 4B as a bounded display preview, not a new biomass grid.
global_file <- file.path(root, "global_l4b_complete.tif")
if (!file.exists(global_file)) {
  ref <- gedifinder("GEDI04_B", ul_lat = 52, ul_lon = -180, lr_lat = -52,
                   lr_lon = 180, version = "002", source = "gee")
  global <- readLevel4B(ref, metrics = "MU")
  map_download(global, global_file, region = c(-180, 180, -52, 52),
               scale = 20000, crs = "EPSG:4326")
}
log("Global Level 4B display saved")
```

### Retrieve and screen several GEDI years

```r
aoi <- st_read(system.file("extdata", "shenandoah.shp", package = "rGEDI"), quiet = TRUE)
box <- st_bbox(aoi)
source_years <- c(2019L, 2021L, 2022L)
all <- list()
for (year in source_years) {
  dest <- file.path(root, paste0("shenandoah_predictors_", year, ".rds"))
  if (file.exists(dest)) {
    all[[as.character(year)]] <- readRDS(dest)
    next
  }
  ref <- gedifinder("GEDI04_A", ul_lat = box[["ymax"]], ul_lon = box[["xmin"]],
    lr_lat = box[["ymin"]], lr_lon = box[["xmax"]], version = "002.1", source = "gee",
    daterange = c(sprintf("%s-05-01", year), sprintf("%s-10-01", year)))
  x <- readLevel4A(ref)
  x <- getLevel4A(x, quality = TRUE, collect = FALSE, max_granules = 500L)
  x <- clip(x, aoi)
  log(year, "collecting park-clipped footprints")
  d <- getLevel4A(x, quality = TRUE, max_features = 10000L, max_granules = 500L,
    cols = c("beam", "shot_number", "agbd", "agbd_se", "sensitivity",
             "solar_elevation", "l4_quality_flag", "degrade_flag"))
  saveRDS(d, file.path(root, paste0("shenandoah_raw_", year, ".rds")))
  d <- clip(d, aoi)
  d <- d[is.finite(agbd) & agbd >= 0 & sensitivity >= .9 & sensitivity <= 1 &
         solar_elevation < 0 & as.character(beam) %in% c("5","6","8","11","BEAM0101","BEAM0110","BEAM1000","BEAM1011")]
  d[, `:=`(sample_id = seq_len(.N) + year*100000L, year = year)]
  log(year, "screened footprints:", nrow(d))
  if (!nrow(d)) next
  predictors <- ee_build_AlphaEarth_embedding_terrain_stack(aoi, year, year, add_lonlat = FALSE)
  bands <- c(sprintf("A%02d", 0:63), "elevation", "slope", "aspect")
  coordinates <- d[,.(sample_id,lon_lowestmode,lat_lowestmode)]
  d <- extractEE(predictors$select(as.list(bands)), d, scale = 30, chunk_size = 1000L)
  d <- merge(d,coordinates,by="sample_id",all.x=TRUE,sort=FALSE)
  d <- d[complete.cases(d[, c("agbd", bands), with = FALSE])]
  d[, `:=`(aspect_sin = sin(aspect*pi/180), aspect_cos = cos(aspect*pi/180))]
  saveRDS(d, dest)
  all[[as.character(year)]] <- d
  log(year, "predictor records:", nrow(d))
}
```

### Select spatial samples and validation blocks

```r
set.seed(726)
all <- lapply(all, function(d) rGEDI::sample(d, spacedSampling(size = min(1000L, nrow(d)), radius = 60)))
data <- rbindlist(all, fill = TRUE)
stopifnot(nrow(data) >= 150L)
xy <- st_coordinates(st_transform(st_as_sf(data, coords = c("lon_lowestmode", "lat_lowestmode"), crs = 4326), 32617))
data[, block := paste(floor(xy[,1]/2000), floor(xy[,2]/2000), sep = "_")]
set.seed(726)
blocks <- base::sample(unique(data$block), ceiling(.3*length(unique(data$block))))
train <- data[!block %in% blocks]; valid <- data[block %in% blocks]
bands <- c(sprintf("A%02d", 0:63), "elevation", "slope", "aspect_sin", "aspect_cos")
selection_file <- file.path(root, "shenandoah_selection.rds")
if (file.exists(selection_file)) selection <- readRDS(selection_file) else {
  log("Selecting variables using training blocks only")
  selection <- varSel(train[, ..bands], train$agbd, method = "rfe", threshold = 0,
                      seed = 726, ntree = 150, nodesize = 5)
  saveRDS(selection, selection_file)
}
keep <- selection$selvars
fwrite(selection$importance, "readme/results/shenandoah_importance.csv")
fwrite(selection$test, "readme/results/shenandoah_rfe.csv")
log("Training:",nrow(train),"validation:",nrow(valid),"selected predictors:",length(keep))
```

### Train, predict and evaluate the GEE forest

```r
# Add circular aspect predictors to the annual GEE stack.
stack <- ee_build_AlphaEarth_embedding_terrain_stack(aoi, 2021, 2021, add_lonlat = FALSE)
aspect <- stack$select("aspect")$multiply(pi/180)
stack <- stack$addBands(aspect$sin()$rename("aspect_sin"))$addBands(aspect$cos()$rename("aspect_cos"))
# Validate on year-matched predictor values in feature tables. The same server
# forest and selected variables are subsequently used for the map.
small <- function(d) d[, c("lon_lowestmode","lat_lowestmode","agbd","sample_id",keep), with = FALSE]
model <- build_ee_forest(vect_as_ee(to_vect(small(train))), "agbd", keep, trees = 500L, seed = 726L)
classify <- function(d) {
  fc <- vect_as_ee(to_vect(small(d)))$classify(model, "predicted")
  out <- reticulate::py_to_r(fc$select(c("sample_id","agbd","predicted"))$getInfo())
  rbindlist(lapply(out$features, `[[`, "properties"))
}
pred_train <- classify(train); pred_train[, stage := "Training"]
pred_valid <- classify(valid); pred_valid[, stage := "Validation"]
predictions <- rbindlist(list(pred_train, pred_valid))
fwrite(predictions, "readme/results/shenandoah_predictions.csv")
fwrite(rbindlist(lapply(split(predictions, predictions$stage), function(d) {
  x <- as.data.table(fit_metrics(d$agbd, d$predicted, unit = "Mg/ha"))
  x[, stage := d$stage[1]]; x
})), "readme/results/shenandoah_accuracy.csv")
fwrite(data[, .(sample_id, year, beam, shot_number, lon_lowestmode, lat_lowestmode, agbd, block,
                 partition = fifelse(block %in% blocks, "Validation", "Training"))],
       "readme/results/shenandoah_footprints.csv")
saveRDS(list(training = data, selected = keep, heldout_blocks = blocks), file.path(root,"shenandoah_model_inputs.rds"))
```

### Refit and export a wall-to-wall map

```r
# Refit on all available observations for the final 2021 demonstration map.
final <- build_ee_forest(vect_as_ee(to_vect(small(data))), "agbd", keep, trees = 500L, seed = 726L)
image <- map_create(final, stack$select(as.list(keep)), aoi = aoi, name = "AGBD")
dest <- file.path(root, "shenandoah_agbd_2021.tif")
if (!file.exists(dest)) map_download(image, dest, aoi, scale = 100, crs = "EPSG:32617")
dem_dest <- file.path(root,"shenandoah_dem.tif")
if (!file.exists(dem_dest)) map_download(stack$select("elevation"),dem_dest,aoi,scale=100,crs="EPSG:32617")
log("Wall-to-wall map saved:", dest)
fwrite(data.table(n = nrow(data), train = nrow(train), validation = nrow(valid),
                  predictors = length(keep), map_year = 2021, map_resolution_m = 100),
       "readme/results/shenandoah_run.csv")
log("COMPLETE")
```

### Inspect the 1:1 plots

Points use a viridis palette to show local footprint counts in 25 × 25 Mg/ha bins. Both panels use
the same axes and colour scale, including all observations. Predictive R² is
one minus residual sum of squares divided by observed total sum of squares.

```r
library(ggplot2)
theme_set(theme_minimal(base_size=12)+theme(panel.grid.minor=element_blank(),plot.title=element_text(face="bold")))
pred <- fread("readme/results/shenandoah_predictions.csv")
pred[,stage:=factor(stage,levels=c("Training","Validation"))]
summary <- pred[,.(n=.N,R2=1-sum((predicted-agbd)^2)/sum((agbd-mean(agbd))^2),
  RMSE=sqrt(mean((predicted-agbd)^2)),bias=mean(predicted-agbd)),by=stage]
summary[,label:=sprintf("n = %d\nR\u00b2 = %.2f\nRMSE = %.1f Mg/ha\nBias = %.1f Mg/ha",n,R2,RMSE,bias)]
maximum <- ceiling(max(pred$agbd,pred$predicted)/100)*100
pred[, `:=`(bin_x=floor(agbd/25),bin_y=floor(predicted/25))]
pred[, local_count := .N, by=.(stage,bin_x,bin_y)]
p <- ggplot(pred[order(local_count)],aes(agbd,predicted))+
  geom_point(aes(colour=local_count),size=1.25,alpha=.9)+
  geom_abline(slope=1,intercept=0,linetype=2,colour="#5D666A",linewidth=.7)+
  scale_colour_viridis_c(option="viridis",trans="log1p",name="Local footprint count",
    breaks=c(1,2,5,10,20,50),guide=guide_colourbar(barwidth=grid::unit(7,"cm"),barheight=grid::unit(.4,"cm")))+
  geom_text(data=summary,aes(x=.04*maximum,y=.96*maximum,label=label),inherit.aes=FALSE,hjust=0,vjust=1,size=4)+
  facet_wrap(~stage)+coord_equal(xlim=c(0,maximum),ylim=c(0,maximum),expand=FALSE)+
  labs(x="Reference GEDI AGBD (Mg/ha)",y="Predicted AGBD (Mg/ha)",title="Biomass model training and spatial validation",
    subtitle="Points coloured by local footprint density (25 x 25 Mg/ha bins)",
    caption="Dashed line: 1:1. Statistics use all observations. Validation holds out 2 km blocks across all sampled years.")+
  theme(strip.text=element_text(face="bold",size=13),panel.background=element_rect(fill="#FAFAF9",colour=NA),
    panel.grid.major=element_line(colour="#E5E7E5"),legend.position="bottom",plot.caption=element_text(hjust=0,size=9))
ggplot2::ggsave(file.path(outdir, "model-validation.png"), p, width = 10.5, height = 6.3, dpi = 300, bg = "white")
output <- summary[,.(stage,n,R2=round(R2,3),RMSE=round(RMSE,2),bias=round(bias,2))]
writeLines(capture.output(print(output)),file.path(outdir, "model-evaluation.txt"))
p
```

```text
        stage     n    R2  RMSE  bias
       <fctr> <int> <num> <num> <num>
1:   Training  1703 0.799 56.03  0.16
2: Validation   645 0.516 88.47 10.34
```

<p align="center"><img src="readme/fig-gee-model-validation.png" width="1000" alt="Training and spatial validation with a shared viridis footprint-count legend"></p>

The displayed run uses 2,348 samples, with 1,703 for training and 645 for
validation, and selects 13 predictors. Validation measures agreement with
withheld GEDI estimates; it is not an independent field biomass validation.

```r
library(patchwork)
importance <- fread("readme/results/shenandoah_importance.csv")
rfe <- fread("readme/results/shenandoah_rfe.csv")
top <- importance[order(-importance)][1:15]
top[,parameter:=factor(parameter,levels=rev(parameter))]
chosen <- sum(importance$selected)
a <- ggplot(top,aes(importance,parameter,fill=selected))+geom_col(width=.7)+
  scale_fill_manual(values=c("FALSE"="#C9D0D0","TRUE"="#705196"),labels=c("Other","Selected"),name=NULL)+
  labs(title="A  Predictor importance",x="Relative importance",y=NULL)+theme(legend.position="bottom")
b <- ggplot(rfe,aes(nvariables,oob_rmse))+geom_line(colour="#267F82",linewidth=.8)+
  geom_vline(xintercept=chosen,colour="#A42A73",linetype=2)+
  labs(title="B  Recursive feature elimination",x="Predictors retained",y="Training OOB RMSE (Mg/ha)")
selection_plot <- a+b+patchwork::plot_annotation(caption="Variable selection uses training blocks only.")
ggplot2::ggsave(file.path(outdir,"variable-selection.png"),selection_plot,width=10.5,height=4.5,dpi=300,bg="white")
selection_plot
```

<p align="center"><img src="readme/fig-gee-variable-selection.png" width="1000" alt="Training variable importance and recursive feature elimination"></p>

### View the final prediction

```r
map_view(list(AGBD = image),
  vis = list(AGBD = list(min = 0, max = 500,
    palette = c("#000004", "#721F81", "#F1605D", "#FCFDBF"))), aoi = aoi)
```

<p align="center"><img src="readme/fig-gee-shenandoah-map.png" width="1050" alt="Shenandoah boundary, GEDI spatial partitions and wall-to-wall biomass map"></p>

The final map refits the selected model on all samples and uses the 2021
predictor stack. It is a prediction between footprints, separate from the
native NASA Level 4B product.

## Close HDF5 connections

```r
close(level1b)
close(level2a)
close(level2b)
close(level2a_h5_clip)
close(sim_clean)
close(sim_noisy)
```

## Reproduce the displayed outputs

Run scripts from the repository root after installing their plotting
dependencies. The local scripts use bundled data; the cloud script is optional.

```r
# install.packages(c("lidR", "plot3D", "patchwork", "viridisLite", "rnaturalearthdata"))
source("readme/build-tutorial-local.R")    # Waveforms, measured ALS and R output
source("readme/build-package-spatial.R")   # Footprint clipping and user grids
source("readme/build-tutorial-rasters.R")  # Native Level 3/4B plots and output
# source("readme/build-gee-tutorial.R")    # Requires GEE; downloads and models
# source("readme/plot-gee-tutorial.R")     # Draw the saved cloud results
```

All figures and data here are package tutorial examples. Data-source and
subset details are recorded in `inst/extdata/EXAMPLE_DATA.txt`.

## Acknowledgements and funding

GEDI data are provided by NASA's LP DAAC and ORNL DAAC. Steven Hancock
developed the GEDI waveform simulator; rGEDI builds on his
[reference implementation](https://bitbucket.org/StevenHancock/gedisimulator/src/master/).
Example boundaries are from the [USDA Forest Service](https://apps.fs.usda.gov/arcx/rest/services/EDW/EDW_ForestSystemBoundaries_01/MapServer/0)
and the [National Park Service](https://www.nps.gov/subjects/gisandmapping/npmap.htm).

This work was primarily supported by NASA’s ICESat-2 Science Team
(80NSSC23K0941), Carbon Monitoring System (80NSSC23K1257), and Commercial
Smallsat Data Scientific Analysis program (80NSSC24K0055); the
McIntire–Stennis Program at the University of Florida (Accession No. 7005758);
and the Joint Fire Science Program through the EMS4D project (22-2-02-15).

## Cite rGEDI

```r
citation("rGEDI")
```

Silva, C. A., and Hamamura, C. *rGEDI: NASA's Global Ecosystem Dynamics
Investigation (GEDI) Data Visualization and Processing*. R package.
[CRAN package DOI](https://doi.org/10.32614/CRAN.package.rGEDI).

## Questions and issues

Report reproducible problems at the
[issue tracker](https://github.com/carlos-alberto-silva/rGEDI/issues), including
the package version, operating system, product/version and minimal example.
Keep credentials out of reports. rGEDI is an independent open-source project,
not an official NASA software product.

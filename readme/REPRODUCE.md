# Reproduce the repository figures

Run these commands from the rGEDI repository root. GEDI access, extraction,
clipping, summary statistics and modelling use rGEDI. Plotting uses R graphics
packages or Python Matplotlib. The displayed observations are measured data.

## Requirements

Install the development rGEDI package and the R packages data.table, sf, terra,
ggplot2, patchwork and MASS. The Python plotting scripts require numpy,
matplotlib and Pillow. Python is needed only to redraw the function diagram
and the combined ALS and waveform plate. Earth Engine access is optional
for the saved examples.

## Function diagram

```text
python readme/build-function-workflow.py
```

This draws the original function diagram as PNG and vector PDF. Labels with
slashes abbreviate related functions. The package reference manual contains
their full names and arguments.

## Forest point clouds and real waveforms

```text
python readme/build-us-forest-figure.py
```

The inputs in `readme/data/us_forests/` contain measured ALS points, DTM
samples, archived Level 1B waveform bins and site metadata. The supplied
`fig-us-forest-sites.png` is the locator map. Rebuild it with
`Rscript readme/build-us-locator.R` from the archived state and Level I
ecoregion vectors. These vectors are transformed to NAD83 Conus Albers
(EPSG 5070) for plotting. Rebuild the locator before the combined plate. Each site has its original
GEDI granule, shot identifier, acquisition date, coordinates and ALS year.
ALS points use height above local ALS ground for colour. Their 3D elevation
is relative to the local reference. Waveforms use height above the independent
GEDI ground estimate. RH98 is the Level 2A value for the recorded shot.

The DTM hillshade uses azimuth 315 degrees, illumination elevation 40 degrees
and no vertical exaggeration. The translucent DTM and 35.264 degree orthographic viewing
elevation keep points visible near the ground. All six clouds share a 30 by 30 m horizontal
view, a -10 to 60 m vertical range, camera angle and 0 to 60 m height colours.
Each waveform uses the same -10 to 60 m height and 0 to 100 percent peak
return scales. Percent peak amplitude is not cumulative waveform energy. The top axis
shows a trapezoidal integration of the nonnegative return over the common
displayed -10 to 60 m interval, normalized to 0 to 1. Official RH98 remains
the L2A value. The cumulative curve does not recompute the mission RH metrics.
The displayed ALS points and DTM are clipped to a 12.5 m radius, with a
translucent surface behind the returns.

The ALS projects are WA Olympic Peninsula 2013, AZ North Kaibab NF B1 2019,
ME Crown of Maine B2 2018, TN Eastern 2 16 B16 Del2 2016, CA Sierra Nevada 6
2022 and FL 2017 Palm Coast 1 C22. The regional labels describe contrasting
forest landscapes. They do not constitute a complete U.S. biome sample. ALS and GEDI
acquisitions differ in year at five sites. These are structural illustrations,
not temporally matched validation measurements.

## Clipping and vegetation profiles

```r
source("readme/build-profile-example.R")
```

The 145 kB `readme/data/appalachian_level2b.h5` is a local subset of NASA
GEDI02_B_2019151032858_O02629_03_T03570_02_004_01_V003, beam BEAM0101,
acquired on 31 May 2019. Its original granule name is recorded in its HDF5
attributes. The script reads it with readLevel2B() and extracts its native
5 m PAVD bins with getLevel2BPAVDProfile(). It screens Level 2B quality and
algorithm flags, degradation, sensitivity and valid cover and PAI values.

Panel A contains 43,638 screened footprints from 13 GEE L2A V2 granules
between 1 May and 1 October 2019. The named NPS park boundary identifies
which footprints are retained by clip(). Panel B contains 4,387 geolocated
footprints across all eight beams of one Earthdata L2A V3 granule, orbit
02629 part 03. Its geometry is extracted with getTrack(). Showing all eight
tracks is a geolocation demonstration, including coverage and power beams.

The selected L2B V3 beam has 686 candidates and 167 screened observations.
Of these, 166 lie inside the park and match valid L2A V3 RH98 values by shot
ID. Their forest classes are sampled from NLCD 2019 with rGEDI extractEE().
There are 72 deciduous, 18 evergreen, 57 mixed-forest and 19 other-cover shots.
The forest comparisons omit the other-cover group. Their means and
interquartile bands use the native PAVD bins and include zeros.

The longest continuous screened section has 43 shots spanning 2.39 km.
PAVD columns occupy display cells ending halfway between adjacent footprint
centres. Vertical interpolation softens each column. No new observations are
generated, and native bins are retained for every summary. L2B ground and
matched L2A RH98 elevation are overlaid. The reserved margin contains a real
925-bin L1B waveform from shot 26290500300273355, marked in B and C.

The boundary is archived in `readme/data/grsm_boundary.geojson`. It comes from
the NPS GRSM_PARK_BOUNDARY service, BOUNDARY_PY layer, with 10 m display
simplification. The original service is recorded in refresh-profile-context.R.
The data sources are the [NPS GIS page](https://www.nps.gov/grsm/planyourvisit/access_data.htm)
and the [USGS NLCD 2019 catalog](https://developers.google.com/earth-engine/datasets/catalog/USGS_NLCD_RELEASES_2019_REL_NLCD).
No background imagery is used in the characterization figure.

All summaries describe sampled footprints, not design-based estimates of
park-wide means. The multigranule overview is GEE V2. Numeric summaries,
the eight-track granule and the waveform/profile pair use Earthdata V3.
Product releases are kept explicit rather than treated as interchangeable.

## Refreshing measured inputs

Cached CSV and vector files support offline redrawing. To retrieve the new
forest examples afresh, configure Earthdata and GEE with rGEDI, set
EE_PROJECT and RGEDI_FIGURE_DATA to a large-data directory, and run
`readme/collect-forest-examples.R`. GEDI selection and extraction use rGEDI.
The measured 3DEP nodes are decoded by collect-als-examples.py and
als_3dep_helpers.py, which require laspy, lazrs, requests, numpy and pyproj.
Use collect-track-overview.R to rebuild the multigranule and eight-track
inputs and the waveform inset. refresh-profile-context.R retrieves the NPS
boundary and NLCD classes. Large retrieval caches default to the D drive.

## Forest biomass model

```r
Sys.setenv(RGEDI_FOREST_OUTPUT = "D:/rGEDI_example/forest_model")
source("inst/scripts/forest_biomass_workflow.R")
source("readme/build-forest-example.R")
```

The default fixture `inst/extdata/forest_model_training.rds` preserves the
original numeric precision and recorded spatial fold assignments. The CSV
version permits inspection. Refitting the rounded CSV can change random
forest splits. Use the RDS to reproduce the illustrated fit.

The dataset pools 392 GEDI L4A V2.1 biomass samples from 2019, 499 from 2021
and 486 from 2022. Each sample has year-matched AlphaEarth predictors and
SRTM terrain covariates extracted through rGEDI. Spatial folds use 2.5 km
blocks. Each validation fit trains on 2019/2021 observations outside its
withheld blocks and predicts 2022 observations inside those blocks. The
pooled holdout has 486 observations. A separate final model uses all 1,377
samples. The script fixes the random seeds and fits 500 trees per model.

Validation gives predictive R2 0.55885, RMSE 56.45299 Mg/ha and bias
4.57947 Mg/ha. The response is GEDI AGBD, so these statistics assess agreement
with held-out GEDI estimates. Both scatterplot panels share axes, density
colours, a fitted line and a 1 to 1 line. All observations contribute to the
statistics. One training observation exceeds the 0 to 750 Mg/ha display,
and a separate full-range plot shows every observation.

## Forest maps and imagery

Set `RGEDI_FOREST_MAP=TRUE` and `EE_PROJECT` before running the forest
workflow to retrieve optional GEE predictors and Landsat imagery through
rGEDI. The direct download example uses 100 m pixels to keep the request
small. The displayed saved biomass map was generated from the same final
reference model on 25 m predictor tiles.

```text
Rscript readme/plot-forest-maps.R DATA_DIR BOUNDARY_GEOJSON
```

DATA_DIR must contain `agbd_2019.tif` and `landsat_2019.tif` in the same
projected metric CRS. The imagery has red, green and blue named bands.
The map helper uses the common biomass extent for both views and adds the
site boundary, north arrow, scale and legends. Biomass units are Mg/ha.
Landsat natural colour uses a fixed surface-reflectance stretch.

## Editable PowerPoint figures

The four files in `readme/editable/` contain the function workflow, the six
forest examples, the linked footprint/profile example and the biomass model
and maps. Workflow boxes, arrows and labels are native editable shapes.
Waveform paths, cumulative curves, RH98 lines, ground references and labels
are editable. The model scatterplots are native charts with embedded data.
Dense point-cloud and PAVD panels are separate image objects. Their internal
plot formatting can be regenerated with the supplied R and Python scripts.

The optional build-editable-figures.mjs uses the Codex bundled presentation
runtime. It does not add a dependency to rGEDI. The canonical scientific PNG
and PDF figures are rebuilt with the R and Python plotting scripts above.
prepare-editable-figure-data.py MODEL_OUTPUT_DIR updates the exact waveform
and model values cached for the deck builder. A new PowerPoint version must
use a new output filename. Existing final files are preserved.

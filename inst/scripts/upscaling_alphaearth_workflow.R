# End-to-end GEDI Level 4A biomass upscaling with AlphaEarth and terrain.
#
# Authentication is read from the user's standard configuration. This script
# never contains or writes NASA Earthdata or Google credentials.

repos <- c(
  rgedi = "https://carlos-alberto-silva.r-universe.dev",
  CRAN = "https://cloud.r-project.org"
)
required <- c("rGEDI", "data.table", "sf", "terra", "reticulate",
              "randomForest", "leaflet")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = repos, dependencies = TRUE)

suppressPackageStartupMessages({
  library(rGEDI)
  library(data.table)
  library(sf)
  library(terra)
})

outdir <- file.path(tempdir(), "rGEDI-alphaearth")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Study area and model period -------------------------------------------------
aoi_path <- system.file("extdata", "stands_cerrado.shp", package = "rGEDI")
aoi <- sf::st_read(aoi_path, quiet = TRUE)
aoi <- sf::st_make_valid(sf::st_transform(aoi, 4326))
box <- sf::st_bbox(aoi)
study_extent <- unname(box[c("xmin", "xmax", "ymin", "ymax")])
daterange <- c("2019-04-18", "2019-04-19")
start_year <- 2019
end_year <- 2019

# Package and service configuration ------------------------------------------
rGEDI_configure(install = TRUE)
earthdata_login() # reads NETRC or ~/.netrc; credentials are never stored here
ee_project <- Sys.getenv("EE_PROJECT", unset = "")
if (!nzchar(ee_project)) {
  stop("Set EE_PROJECT to your Google Cloud project ID before running.")
}
ee <- ee_initialize(project = ee_project)

# Search and stream Level 4A --------------------------------------------------
urls <- gedifinder(
  "GEDI04_A",
  ul_lat = box[["ymax"]], ul_lon = box[["xmin"]],
  lr_lat = box[["ymin"]], lr_lon = box[["xmax"]],
  version = "003", daterange = daterange,
  cloud_computing = TRUE, persist = TRUE
)
level4a <- readLevel4A(urls[[1L]])
on.exit(close(level4a), add = TRUE)

footprints <- getLevel4A(
  level4a,
  cols = c(
    "beam", "shot_number", "delta_time", "lat_lowestmode",
    "lon_lowestmode", "agbd", "agbd_se", "sensitivity",
    "l4a_quality_flag_rel3", "degrade_include_flag",
    "elev_highestreturn_outlier_flag"
  ),
  quality = TRUE
)
footprints <- clipLevel4AGeometry(footprints, aoi)
footprints <- footprints[is.finite(agbd) & agbd >= 0]

# Spatially balanced sampling and predictor extraction -----------------------
set.seed(42)
sampled <- sampleGEDI(
  footprints,
  spacedSampling(size = min(500L, nrow(footprints)), radius = 25)
)
sampled$sample_id <- seq_len(nrow(sampled))

stack <- ee_build_AlphaEarth_embedding_terrain_stack(
  aoi, start_year = start_year, end_year = end_year, add_lonlat = TRUE
)
bands <- reticulate::py_to_r(stack$bandNames()$getInfo())
predictor_names <- c(
  grep("^A[0-9]{2}$", bands, value = TRUE),
  intersect(c("elevation", "slope", "aspect", "longitude", "latitude"), bands)
)

training <- extractEE(
  stack$select(as.list(predictor_names)), sampled,
  scale = 30, chunk_size = 250
)
coordinates <- as.data.table(sampled)[, .(
  sample_id, lon_lowestmode, lat_lowestmode
)]
if (!all(c("lon_lowestmode", "lat_lowestmode") %in% names(training))) {
  training <- merge(training, coordinates, by = "sample_id", all.x = TRUE,
                    sort = FALSE)
}
complete <- training[
  complete.cases(training[, c("agbd", predictor_names), with = FALSE])
]
fwrite(complete, file.path(outdir, "gedi-alphaearth-training.csv"))

# False-color AlphaEarth visualization ---------------------------------------
rgb <- stack$select(c("A00", "A20", "A40"))
rgb_map <- map_view(
  list(`AlphaEarth RGB` = rgb),
  vis = list(`AlphaEarth RGB` = list(
    bands = c("A00", "A20", "A40"), min = -0.06, max = 0.12
  )),
  aoi = aoi
)
print(rgb_map)

# RFE and independent 30% holdout --------------------------------------------
selection <- varSel(
  complete[, ..predictor_names], complete$agbd,
  method = "rfe", threshold = 0, seed = 42, ntree = 200
)
best_predictors <- selection$selvars
plot(selection, which = "importance", main = "RFE predictor importance")
plot(selection, which = "rfe", main = "RFE out-of-bag error")

fit <- fit_model(
  complete[, ..best_predictors], complete$agbd,
  method = "randomForest", test = "split", test_size = 0.30,
  seed = 42, ntree = 500, importance = TRUE
)
print(fit$stats_train)
print(fit$stats_test)

# Train the equivalent Earth Engine regressor and map the AOI ----------------
ee_columns <- c("agbd", "lon_lowestmode", "lat_lowestmode", best_predictors)
training_vect <- to_vect(
  complete[, ..ee_columns], lon = "lon_lowestmode", lat = "lat_lowestmode"
)
training_ee <- vect_as_ee(training_vect)$filter(
  ee$Filter$notNull(as.list(c("agbd", best_predictors)))
)
forest <- build_ee_forest(
  training_ee, response = "agbd", predictors = best_predictors,
  trees = 500, seed = 42
)
agbd_map <- map_create(
  forest, stack$select(as.list(best_predictors)),
  aoi = aoi, name = "agbd"
)

agbd_view <- map_view(
  list(`Predicted AGBD` = agbd_map),
  vis = list(`Predicted AGBD` = list(
    min = 0, max = 200,
    palette = c("#f7fcf5", "#74c476", "#00441b")
  )),
  aoi = aoi
)
print(agbd_view)

# Export large maps through Earth Engine/Google Drive -------------------------
drive_task <- ee_image_to_drive(
  agbd_map,
  description = "rGEDI_AGBD_2019",
  folder = "EE_Exports",
  file_name_prefix = "rGEDI_AGBD_2019",
  region = aoi,
  scale = 30,
  start = TRUE
)
print(ee_check_task_status(drive_task, quiet = FALSE))

# For an AOI small enough for a direct request:
# map_download(
#   agbd_map, file.path(outdir, "rGEDI_AGBD_2019.tif"),
#   region = aoi, scale = 30, overwrite = TRUE
# )

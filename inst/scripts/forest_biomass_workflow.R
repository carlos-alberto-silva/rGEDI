# Forest modelling example for rGEDI 0.6.1.
# Run from the rGEDI repository root, or set RGEDI_FOREST_TRAINING to a saved RDS or CSV table.
# Training responses are GEDI L4A V2.1 estimates for 2019/2021/2022. Annual
# AlphaEarth and SRTM predictors were extracted with rGEDI's GEE interface.
library(rGEDI)
library(data.table)
library(sf)
library(terra)
input <- Sys.getenv("RGEDI_FOREST_TRAINING",
                   unset = system.file("extdata", "forest_model_training.rds", package = "rGEDI"))
if (!nzchar(input) || !file.exists(input)) {
  input <- "inst/extdata/forest_model_training.rds"
}
outdir <- Sys.getenv("RGEDI_FOREST_OUTPUT", unset = file.path(tempdir(), "rGEDI-forest"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
d <- if (grepl("[.]rds$", input, ignore.case = TRUE)) readRDS(input) else fread(input)
predictor_names <- c(sprintf("A%02d", 0:63), "srtm_elevation", "slope", "aspect_sin", "aspect_cos")
stopifnot(all(c("agbd", "year", "lon", "lat", "fold", predictor_names) %in% names(d)))
training <- copy(d)
validation <- list()
for (f in sort(unique(d$fold))) {
  train <- d[fold != f & year < 2022]
  holdout <- d[fold == f & year == 2022]
  if (!nrow(holdout)) next
  fit <- fit_model(train[, ..predictor_names], train$agbd,
                   method = "randomForest", test = "none", seed = 5721 + f,
                   ntree = 500, mtry = floor(sqrt(length(predictor_names))), nodesize = 5)
  validation[[length(validation) + 1L]] <- data.table(
    observed = holdout$agbd,
    predicted = as.numeric(predict(fit, newdata = holdout[, ..predictor_names])),
    fold = f, stage = "Spatial and 2022 holdout")
}
model <- fit_model(d[, ..predictor_names], d$agbd, method = "randomForest", test = "none",
                   seed = 5721, ntree = 500,
                   mtry = floor(sqrt(length(predictor_names))), nodesize = 5)
results <- rbindlist(c(list(data.table(observed = d$agbd, predicted = as.numeric(model$fitted),
                                      fold = NA_integer_, stage = "Training")), validation))
fwrite(results, file.path(outdir, "predictions.csv"))
saveRDS(model, file.path(outdir, "agbd_model.rds"))
print(results[, .(n = .N, R2 = 1 - sum((observed - predicted)^2) / sum((observed - mean(observed))^2),
                  RMSE = sqrt(mean((observed - predicted)^2)), bias = mean(predicted - observed)), by = stage])
print(fit_metrics(results[stage == "Spatial and 2022 holdout"]$observed,
                  results[stage == "Spatial and 2022 holdout"]$predicted, unit = "Mg/ha"))

# Optional GEE predictor download and map prediction.
# Set RGEDI_FOREST_MAP=TRUE and EE_PROJECT to enable it.
if (identical(Sys.getenv("RGEDI_FOREST_MAP"), "TRUE")) {
  project_id <- Sys.getenv("EE_PROJECT")
  if (!nzchar(project_id)) stop("Set EE_PROJECT before requesting the map inputs.")
  ee <- ee_initialize(project = project_id, authenticate = FALSE)
  study_extent <- c(-112.48, -112.20, 36.54, 36.77)
  stack <- ee_build_AlphaEarth_embedding_terrain_stack(study_extent, 2019, 2019, add_lonlat = FALSE)
  path <- file.path(outdir, "predictors_2019.tif")
  if (!file.exists(path)) {
    map_download(stack, path, region = study_extent, scale = 100, crs = "EPSG:32612")
  }
  predictors <- rast(path)
  names(predictors) <- c(sprintf("A%02d", 0:63), "srtm_elevation", "slope", "aspect")
  predictors$aspect_sin <- sin(predictors$aspect * pi / 180)
  predictors$aspect_cos <- cos(predictors$aspect * pi / 180)
  mapped <- terra::predict(predictors[[predictor_names]], model, na.rm = TRUE,
                            filename = file.path(outdir, "agbd_2019.tif"), overwrite = TRUE)
  reflectance <- function(image) {
    valid <- image$select("QA_PIXEL")$bitwiseAnd(63L)$eq(0L)$And(image$select("QA_RADSAT")$eq(0L))
    image$select(list("SR_B4", "SR_B3", "SR_B2"), list("red", "green", "blue"))$
      multiply(.0000275)$add(-.2)$updateMask(valid)
  }
  landsat <- ee$ImageCollection("LANDSAT/LC08/C02/T1_L2")$
    filterBounds(ext_to_ee(study_extent))$filterDate("2019-08-01", "2019-10-01")$
    map(reflectance)$median()
  map_download(landsat, file.path(outdir, "landsat_2019.tif"),
                region = study_extent, scale = 100, crs = "EPSG:32612", overwrite = TRUE)
  event <- ee$FeatureCollection("USFS/GTAC/MTBS/burned_area_boundaries/v1")$
    filter(ee$Filter$eq("Event_ID", "AZ3662811235020200608"))
  info <- reticulate::py_to_r(event$getInfo())
  jsonlite::write_json(info, file.path(outdir, "boundary.geojson"), auto_unbox = TRUE, digits = NA)
}

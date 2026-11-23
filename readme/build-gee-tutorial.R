# Optional package tutorial. Requires an authenticated Earth Engine project.
# Large intermediate data stay in RGEDI_EXAMPLE_CACHE, outside the package.
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(rGEDI)
  library(sf); library(terra); library(data.table)
})
main <- function() {
  root <- Sys.getenv("RGEDI_EXAMPLE_CACHE", file.path(tempdir(), "rGEDI-package-examples"))
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

  boundary_file <- file.path(root, "shenandoah_boundary.geojson")
  if (!file.exists(boundary_file)) {
    url <- paste0("https://services1.arcgis.com/fBc8EJBxQRMcHlei/arcgis/rest/services/",
      "NPS_Land_Resources_Division_Boundary_and_Tract_Data_Service/FeatureServer/2/query?",
      "where=UNIT_CODE%20%3D%20%27SHEN%27&outFields=UNIT_CODE,UNIT_NAME&outSR=4326&f=geojson")
    curl::curl_download(url, boundary_file)
  }
  aoi <- st_read(boundary_file, quiet = TRUE) |> st_make_valid() |> st_transform(32617)
  aoi <- st_simplify(aoi, dTolerance = 20, preserveTopology = TRUE) |> st_transform(4326)
  aoi <- aoi[, c("UNIT_CODE", "UNIT_NAME")]
  st_write(aoi, "inst/extdata/shenandoah.shp", delete_layer = TRUE, quiet = TRUE)
  st_write(aoi, "readme/shenandoah-boundary.geojson", delete_dsn = TRUE, quiet = TRUE)
  box <- st_bbox(aoi)
  log("Boundary saved; area km2:", round(sum(st_area(st_transform(aoi, 32617)))/1e6, 1))
  source_years <- c(2019L, 2021L, 2022L)
  all <- list()
  for (year in source_years) {
    dest <- file.path(root, paste0("shenandoah_predictors_", year, ".rds"))
    if (file.exists(dest)) {
      cached <- readRDS(dest)
      if (!"lon_lowestmode" %in% names(cached)) {
        raw <- readRDS(file.path(root,paste0("shenandoah_raw_",year,".rds")))
        i <- match(cached$shot_number,raw$shot_number)
        stopifnot(!anyNA(i))
        cached[, `:=`(lon_lowestmode=raw$lon_lowestmode[i],lat_lowestmode=raw$lat_lowestmode[i])]
        saveRDS(cached,dest)
      }
      all[[as.character(year)]] <- cached; next
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
}
main()

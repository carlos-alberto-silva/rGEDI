# Rebuild the real-data figures and animation used by README.md.
#
# Requirements:
#   1. NASA Earthdata credentials in NETRC or a standard ~/.netrc file.
#   2. Earth Engine access to the project in EE_PROJECT.
#   3. Optional Python modules installed with rGEDI_configure(install = TRUE).

if (file.exists("DESCRIPTION") && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(rGEDI)
}

ee_project <- Sys.getenv("EE_PROJECT", unset = "")
if (!nzchar(ee_project)) {
  stop("Set EE_PROJECT to your Google Cloud project before running this script.")
}
earthdata_login() # reads NETRC or ~/.netrc; credentials are never stored here

xmin <- -44.18
xmax <- -44.05
ymin <- -13.78
ymax <- -13.67
daterange <- c("2019-04-18", "2019-04-19")
study_extent <- c(xmin, xmax, ymin, ymax)

urls <- gedifinder(
  "GEDI04_A",
  ul_lat = ymax, ul_lon = xmin, lr_lat = ymin, lr_lon = xmax,
  version = "003", daterange = daterange,
  cloud_computing = TRUE, persist = TRUE
)

level4a <- readLevel4A(urls[1])
biomass_orbit <- getLevel4A(
  level4a,
  cols = c(
    "beam", "shot_number", "delta_time", "lat_lowestmode",
    "lon_lowestmode", "agbd", "agbd_se", "sensitivity",
    "l4a_quality_flag_rel3", "degrade_include_flag",
    "elev_highestreturn_outlier_flag"
  ),
  quality = TRUE,
  beams = "BEAM0000"
)
close(level4a)

biomass <- clipLevel4A(biomass_orbit, xmin, xmax, ymin, ymax)
dir.create("readme", showWarnings = FALSE)

orbit_plot <- biomass_orbit[seq(
  1, nrow(biomass_orbit),
  length.out = min(4000, nrow(biomass_orbit))
)]
plot_gedi_orbit_animation(
  orbit_plot,
  output_file = "readme/gedi-orbit-animation.html",
  title = "GEDI Level 4A V3 biomass orbit - BEAM0000",
  duration = 20,
  launch = FALSE
)

palette <- grDevices::colorRampPalette(
  c("#2c115f", "#1fa187", "#fde725")
)(100)

grDevices::png("readme/fig-gedi-cloud-l4a.png", 1600, 850, res = 150)
graphics::par(
  mfrow = c(1, 2), mar = c(4.2, 4.4, 3, 1.2),
  bg = "#07131d", fg = "white", col.axis = "white",
  col.lab = "white", col.main = "white"
)
graphics::plot(
  orbit_plot$lon_lowestmode, orbit_plot$lat_lowestmode,
  pch = 16, cex = 0.18,
  col = grDevices::adjustcolor("#4dd5a5", 0.65),
  xlab = "Longitude", ylab = "Latitude", main = "Streamed GEDI orbit"
)
graphics::grid(col = "#ffffff25")
z <- pmax(
  0,
  pmin(biomass$agbd, stats::quantile(biomass$agbd, 0.98, na.rm = TRUE))
)
index <- pmax(
  1L,
  pmin(100L, 1L + floor(99 * (z - min(z)) / max(1e-9, diff(range(z)))))
)
graphics::plot(
  biomass$lon_lowestmode, biomass$lat_lowestmode,
  pch = 16, cex = 1.5, col = palette[index],
  xlab = "Longitude", ylab = "Latitude",
  main = "Quality-filtered AGBD footprints"
)
graphics::grid(col = "#ffffff25")
grDevices::dev.off()

agbd_raster <- rasterizeLevel4A(
  biomass, metric = "agbd", res = 0.002, fun = mean
)
grDevices::png("readme/fig-gedi-l4a-raster.png", 1200, 900, res = 150)
terra::plot(
  agbd_raster, col = palette,
  main = "GEDI Level 4A mean AGBD (Mg/ha)",
  axes = TRUE, plg = list(title = "AGBD")
)
graphics::points(
  biomass$lon_lowestmode, biomass$lat_lowestmode,
  pch = 16, cex = 0.22,
  col = grDevices::adjustcolor("black", 0.35)
)
grDevices::dev.off()

set.seed(42)
samples <- sampleGEDI(
  biomass, spacedSampling(size = 100, radius = 25)
)
samples$shot_number <- as.character(samples$shot_number)
samples$sample_id <- seq_len(nrow(samples))

ee <- ee_initialize(ee_project, authenticate = FALSE, quiet = TRUE)
full_stack <- ee_build_AlphaEarth_embedding_terrain_stack(
  study_extent, start_year = 2019, end_year = 2019
)
available_bands <- reticulate::py_to_r(full_stack$bandNames()$getInfo())
predictor_names <- c(
  head(grep("^A", available_bands, value = TRUE), 8),
  intersect(c("elevation", "slope", "aspect"), available_bands)
)
stack <- full_stack$select(as.list(predictor_names))

rgb_image <- full_stack$select(c("A00", "A20", "A40"))$visualize(
  bands = c("A00", "A20", "A40"), min = -0.06, max = 0.12
)
map_download(
  rgb_image, "readme/alphaearth-rgb.tif",
  region = study_extent, scale = 30, overwrite = TRUE
)
rgb_raster <- terra::rast("readme/alphaearth-rgb.tif")
grDevices::png("readme/fig-alphaearth-rgb.png", 1200, 900, res = 150)
terra::plotRGB(rgb_raster, r = 1, g = 2, b = 3, stretch = "lin",
               main = "AlphaEarth embedding false-color composite")
grDevices::dev.off()

training <- extractEE(stack, samples, scale = 30, chunk_size = 50)
training <- training[stats::complete.cases(
  training[, c("agbd", predictor_names), with = FALSE]
)]
if (!all(c("lon_lowestmode", "lat_lowestmode") %in% names(training))) {
  coordinates <- data.table::as.data.table(samples)[, .(
    sample_id, lon_lowestmode, lat_lowestmode
  )]
  training <- merge(training, coordinates, by = "sample_id", sort = FALSE)
}
data.table::fwrite(training, "readme/gedi-alphaearth-training.csv")

selection <- varSel(
  training[, predictor_names, with = FALSE], training$agbd,
  method = "rfe", threshold = 0, seed = 42, ntree = 200
)
predictor_names <- selection$selvars
grDevices::png("readme/fig-gedi-rfe.png", 1500, 750, res = 150)
graphics::par(mfrow = c(1, 2), mar = c(4.4, 7.5, 3, 1))
plot(selection, which = "importance", main = "Predictor importance")
plot(selection, which = "rfe", main = "Recursive feature elimination")
grDevices::dev.off()

model <- fit_model(
  training[, predictor_names, with = FALSE], training$agbd,
  method = "randomForest", test = "kfold", k = 5,
  seed = 42, ntree = 100
)

observed <- model$response
predicted <- model$validation
ok <- is.finite(predicted)
grDevices::png(
  "readme/fig-gedi-model-validation.png", 1100, 900, res = 150
)
graphics::par(mar = c(4.5, 4.7, 3.2, 1.2))
graphics::plot(
  observed[ok], predicted[ok], pch = 21,
  bg = "#1fa187", col = "#173f5f",
  xlab = "Observed GEDI AGBD (Mg/ha)",
  ylab = "Five-fold prediction (Mg/ha)",
  main = "Demonstration model validation"
)
graphics::abline(0, 1, lwd = 2, col = "#d1495b")
graphics::grid()
graphics::legend(
  "topleft", bty = "n",
  legend = sprintf("%s = %.3f", model$stats_test$stat,
                   model$stats_test$value)
)
grDevices::dev.off()

training_ee <- vect_as_ee(to_vect(
  training, lon = "lon_lowestmode", lat = "lat_lowestmode"
))
training_ee <- training_ee$filter(
  ee$Filter$notNull(as.list(c("agbd", predictor_names)))
)
ee_model <- build_ee_forest(
  training_ee, "agbd", predictor_names, trees = 100, seed = 42
)
prediction <- map_create(
  ee_model, stack, study_extent, name = "agbd"
)
map_download(
  prediction, "readme/gedi-wall-to-wall-agbd.tif",
  region = study_extent, scale = 30, overwrite = TRUE
)

wall <- terra::rast("readme/gedi-wall-to-wall-agbd.tif")
grDevices::png("readme/fig-gedi-wall-to-wall.png", 1200, 900, res = 150)
terra::plot(
  wall, col = palette,
  main = "GEE wall-to-wall AGBD demonstration",
  axes = TRUE, plg = list(title = "AGBD")
)
graphics::points(
  biomass$lon_lowestmode, biomass$lat_lowestmode,
  pch = 16, cex = 0.28,
  col = grDevices::adjustcolor("white", 0.6)
)
grDevices::dev.off()

print(model$stats_test)
message("README assets rebuilt from ", basename(unclass(urls[1])))

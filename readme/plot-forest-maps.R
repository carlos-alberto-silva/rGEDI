# Draw projected forest rasters with a boundary, north arrow, scale and units.
# Rscript readme/plot-forest-maps.R DATA_DIR BOUNDARY_GEOJSON
library(terra)
library(sf)
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply a directory containing agbd_2019.tif and landsat_2019.tif, plus the boundary GeoJSON.")
folder <- args[1L]
boundary <- sf::st_read(args[2L], quiet = TRUE)
reference <- terra::rast(file.path(folder, "agbd_2019.tif"))
shared_bounds <- terra::ext(reference)
draw <- function(path, optical = FALSE) {
  r <- terra::rast(path)
  if (sf::st_is_longlat(sf::st_crs(terra::crs(r)))) stop("Map inputs must use a projected metric CRS.")
  if (!terra::same.crs(r, reference)) stop("Both map inputs must use the same projected CRS.")
  r <- terra::crop(r, shared_bounds)
  b <- as.vector(shared_bounds)
  fire <- st_transform(boundary, terra::crs(r))
  data <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  if (optical) {
    stopifnot(terra::nlyr(r) >= 3L)
    channels <- as.matrix(data[, c("red", "green", "blue")])
    channels[channels <= -9999] <- NA_real_
    channels[] <- pmin(1, pmax(0, (channels - .015) / .23))^(1 / 1.8)
    keep <- complete.cases(channels)
    data$colour <- NA_character_
    data$colour[keep] <- grDevices::rgb(channels[keep, , drop = FALSE])
    p <- ggplot(data) + geom_raster(aes(x, y, fill = colour)) + scale_fill_identity(na.value = "#EEF0EC")
    title <- "North Kaibab forest landscape"
    key <- "fig-forest-landscape"
    line <- "#F9F3BC"
  } else {
    names(data)[3L] <- "agbd"
    p <- ggplot(data) + geom_raster(aes(x, y, fill = agbd)) +
      scale_fill_viridis_c(option = "magma", name = "AGBD (Mg/ha)", limits = c(0, 250),
                           oob = scales::squish, breaks = seq(0, 250, 50),
                           labels = c("0", "50", "100", "150", "200", "\u2265250"))
    title <- "Modelled aboveground biomass density"
    key <- "fig-gedi-wall-to-wall"
    line <- "#203535"
  }
  # UTM metres are used for the map views and distance labels.
  x <- b[1L] + .07 * diff(b[1:2])
  y <- b[3L] + .07 * diff(b[3:4])
  length_m <- 5000
  p <- p + geom_sf(data = fire, aes(linetype = "2020 Mangum Fire boundary"), fill = NA,
                    colour = line, linewidth = .4) +
    scale_linetype_manual(values = "solid", name = NULL) +
    coord_sf(xlim = b[1:2], ylim = b[3:4], datum = NA, expand = FALSE) +
    annotate("label", x = x + length_m / 2, y = y + .025 * diff(b[3:4]),
             label = "5 km", size = 3, fill = "white", linewidth = 0) +
    annotate("segment", x = x, xend = x + length_m, y = y, yend = y, linewidth = .9) +
    annotate("segment", x = b[2] - .08 * diff(b[1:2]), xend = b[2] - .08 * diff(b[1:2]),
             y = b[4] - .15 * diff(b[3:4]), yend = b[4] - .06 * diff(b[3:4]),
             arrow = grid::arrow(type = "closed", length = grid::unit(2, "mm")), colour = "#203535") +
    annotate("label", x = b[2] - .08 * diff(b[1:2]), y = b[4] - .035 * diff(b[3:4]),
             label = "N", size = 3, fill = "white", linewidth = 0) +
    theme_void(base_size = 10) +
    labs(title = title, subtitle = if (optical) "Landsat 8 natural colour - 2019" else "Mangum forest landscape - 2019") +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
  ggsave(file.path("readme", paste0(key, ".png")), p, width = 7.4, height = 5.4, dpi = 400, bg = "white")
}
draw(file.path(folder, "agbd_2019.tif"))
draw(file.path(folder, "landsat_2019.tif"), optical = TRUE)

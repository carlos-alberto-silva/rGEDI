# Package-only example: original bundled fixtures, no manuscript inputs.
suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".", quiet = TRUE) else library(rGEDI)
  library(data.table); library(sf); library(terra); library(ggplot2); library(patchwork)
})
main <- function() {
  work <- tempfile("rgedi-readme-"); dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  open_product <- function(short, archive) {
    read <- if (short == "2a") readLevel2A else readLevel2B
    read(unzip(system.file("extdata", archive, package = "rGEDI"), exdir = work))
  }
  a <- open_product("2a", "GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.zip")
  b <- open_product("2b", "GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.zip")
  on.exit(close(a), add = TRUE); on.exit(close(b), add = TRUE)
  height <- getLevel2AM(a); vegetation <- getLevel2BVPM(b)
  area <- st_read(system.file("extdata", "stands_cerrado.shp", package = "rGEDI"), quiet = TRUE)
  area$stand_id <- seq_len(nrow(area))
  h <- clip(height, area, split_by = "stand_id")
  v <- clip(vegetation, area, split_by = "stand_id")
  height_summary <- function(x) list(n = length(x), rh98 = mean(x))
  vegetation_summary <- function(pai, cover) list(pai = mean(pai), cover = 100 * mean(cover))
  hs <- polyStatsLevel2AM(h, height_summary(rh98), id = "poly_id")
  vs <- polyStatsLevel2BVPM(v, vegetation_summary(pai, cover), id = "poly_id")
  summary <- merge(hs, vs, by = "poly_id")
  fwrite(summary, "readme/example-stand-summary.csv")
  projected <- st_transform(area, 32723)
  xy <- to_vect(height)
  xy <- terra::project(xy, "EPSG:32723")
  bound <- as.vector(terra::ext(xy))
  bound <- c(floor(bound[1]/1000)*1000, ceiling(bound[2]/1000)*1000,
             floor(bound[3]/1000)*1000, ceiling(bound[4]/1000)*1000)
  template <- terra::rast(terra::ext(bound), resolution = 1000, crs = "EPSG:32723")
  rh <- rGEDI::rasterize(height, "rh98", template = template, fun = "mean")
  pai <- rGEDI::rasterize(vegetation, "pai", template = template, fun = "mean")
  point_sf <- st_as_sf(as.data.frame(height), coords = c("lon_lowestmode", "lat_lowestmode"), crs = 4326) |> st_transform(32723)
  selected <- clip(height, area[1:2, ])
  inside_sf <- st_as_sf(as.data.frame(selected), coords = c("lon_lowestmode", "lat_lowestmode"), crs = 4326) |> st_transform(32723)
  dx <- bound[2]-bound[1]; dy <- bound[4]-bound[3]
  decorate <- function(p) p +
    geom_sf(data = projected, fill = NA, colour = "#53655C", linewidth = .35, inherit.aes = FALSE) +
    coord_sf(xlim = bound[1:2], ylim = bound[3:4], datum = NA, expand = FALSE) +
    annotate("segment", x = bound[1]+.06*dx, xend = bound[1]+.06*dx+2000,
             y = bound[3]+.06*dy, yend = bound[3]+.06*dy, linewidth = .8) +
    annotate("label", x = bound[1]+.06*dx+1000, y = bound[3]+.1*dy, label = "2 km", size = 3.2, linewidth = 0) +
    annotate("segment", x = bound[2]-.09*dx, xend = bound[2]-.09*dx,
             y = bound[4]-.16*dy, yend = bound[4]-.06*dy,
             arrow = grid::arrow(type = "closed", length = grid::unit(2,"mm"))) +
    annotate("text", x = bound[2]-.09*dx, y = bound[4]-.025*dy, label = "N", size = 3.5) +
    theme_void(base_size = 12) + theme(plot.title = element_text(face = "bold"),
      legend.position = "bottom", legend.title = element_text(size = 11),
      plot.background = element_rect(fill = "white", colour = NA))
  first <- ggplot() + geom_sf(data = point_sf, colour = "#CBD1D4", size = .8) +
    geom_sf(data = inside_sf, aes(colour = rh98), size = 1.1) +
    scale_colour_viridis_c(name = "RH98 (m)", limits = range(height$rh98, na.rm = TRUE)) +
    labs(title = "A  Clip footprints", subtitle = "Colour: first two stands; grey: excluded")
  grid_panel <- function(r, title, unit) {
    d <- as.data.frame(r, xy = TRUE, na.rm = TRUE); names(d)[3] <- "value"
    ggplot(d) + geom_tile(aes(x, y, fill = value), width = 1000, height = 1000) +
      scale_fill_viridis_c(name = if (unit == "RH98 (m)") unit else expression(PAI~(m^2/m^2)),
                           limits = if (unit == "RH98 (m)") range(height$rh98, na.rm = TRUE) else NULL) +
      labs(title = title, subtitle = "Mean of footprints in 1 km cells")
  }
  plate <- decorate(first) + decorate(grid_panel(rh, "B  Grid height", "RH98 (m)")) +
    decorate(grid_panel(pai, "C  Grid plant area", "PAI (m2/m2)")) +
    plot_annotation(caption = "Bundled GEDI Level 2A/2B subset, 18 April 2019. Outlines: example stands. White cells: no sampled footprints. UTM 23S.",
      theme = theme(plot.caption = element_text(size = 10, hjust = 0)))
  ggsave("readme/fig-package-spatial.png", plate, width = 12.5, height = 5.4, dpi = 300, bg = "white")
  ggsave("readme/fig-package-spatial.pdf", plate, width = 12.5, height = 5.4, bg = "white")
  print(summary)
}
main()

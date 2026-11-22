# Rebuild the README figures that use rGEDI's bundled GEDI subsets.
library(devtools)
load_all(".", quiet = TRUE)

dir.create("readme", showWarnings = FALSE)
work <- tempfile("rgedi-readme-")
dir.create(work)

level1b_path <- unzip("inst/extdata/GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.zip", exdir = work)
level2a_path <- unzip("inst/extdata/GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.zip", exdir = work)
level2b_path <- unzip("inst/extdata/GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.zip", exdir = work)
stands <- sf::st_read("inst/extdata/stands_cerrado.shp", quiet = TRUE)

level1b <- readLevel1B(level1b_path)
level2a <- readLevel2A(level2a_path)
level2b <- readLevel2B(level2b_path)

geo <- getLevel1BGeo(level1b)
metrics2a <- getLevel2AM(level2a)
metrics2b <- getLevel2BVPM(level2b)
pai <- getLevel2BPAIProfile(level2b)
pavd <- getLevel2BPAVDProfile(level2b)

palette <- grDevices::colorRampPalette(c("#2c115f", "#1fa187", "#fde725"))(100)

# Extracted with getTrack(..., every = 200) from all four Level 2A
# production granules for orbit O01964. The CSV stores coordinates and
# beam/time/granule identifiers only.
track <- if (file.exists("readme/gedi-orbit-track.csv")) {
  data.table::fread("readme/gedi-orbit-track.csv")
} else {
  getTrack(level2a)
}
plot_gedi_orbit_animation(track, output_file = "readme/gedi-orbit-animation.gif",
  title = "GEDI aboard the International Space Station", duration = 8,
  launch = FALSE)
plot_gedi_orbit_animation(track, output_file = "readme/gedi-orbit-animation.html",
  title = "GEDI aboard the International Space Station", duration = 8,
  launch = FALSE)

grDevices::png("readme/fig-gedi-rgt.png", 1200, 850, res = 150)
plot(track$longitude, track$latitude, type = "n",
     xlab = "Longitude", ylab = "Latitude", main = "GEDI reference ground track")
for (beam in unique(track$track)) {
  part <- track[track == beam]
  lines(part$longitude, part$latitude, lwd = 2, col = "#1fa187")
  points(part$longitude, part$latitude, pch = 16, cex = .45, col = "#440154")
}
grid()
grDevices::dev.off()

shot <- metrics2a$shot_number[1]
grDevices::png("readme/fig-waveform-rh.png", 1050, 900, res = 150)
plotWFMetrics(level1b, level2a, shot_number = shot, rh = c(25, 50, 75, 90, 98))
grDevices::dev.off()

row <- which.max(rowSums(as.matrix(pai[, grep("^pai_z", names(pai)), with = FALSE]), na.rm = TRUE))
profile_height <- seq(2.5, 147.5, by = 5)
pai_values <- as.numeric(pai[row, grep("^pai_z", names(pai)), with = FALSE])
pavd_values <- as.numeric(pavd[row, grep("^pavd_z", names(pavd)), with = FALSE])
box <- sf::st_bbox(stands)
clipped2a <- clipLevel2AMGeometry(metrics2a, stands)
clipped2b <- clipLevel2BVPMGeometry(metrics2b, stands)
grid2a <- gridStatsLevel2AM(clipped2a, mean(rh98), res = 0.002)
grid2b <- gridStatsLevel2BVPM(metrics2b, mean(cover), res = 0.002)
set.seed(11)
cloud <- data.frame(
  X = c(runif(1500, -12, 12), runif(3500, -12, 12)),
  Y = c(runif(1500, -12, 12), runif(3500, -12, 12)),
  Z = c(rnorm(1500, 0, .25), pmax(0, rnorm(3500, 17, 6))),
  Classification = c(rep(2L, 1500), rep(5L, 3500))
)
clean <- gediWFSimulator(cloud, output = file.path(work, "sim-clean.h5"),
                         coords = c(0, 0), seed = 11)
noisy <- gediWFSimulator(cloud, output = file.path(work, "sim-noisy.h5"),
                         coords = c(0, 0), noise = 0.03, seed = 11)
data.table::fwrite(metrics2a[, .(beam, shot_number, lon_lowestmode,
  lat_lowestmode, elev_lowestmode, rh50, rh75, rh90, rh98, rh100)],
  "readme/gedi-level2a-example.csv")
close(level1b); close(level2a); close(level2b)
close(clean); close(noisy)
unlink(work, recursive = TRUE)
message("Local README figures and example table rebuilt.")

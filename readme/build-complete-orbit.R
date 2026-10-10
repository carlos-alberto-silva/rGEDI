# Build the README orbit animation from all four GEDI02_A granules that make
# up GEDI orbit O01964. Authentication is read by earthdata_login() from the
# user's normal Earthdata environment or ~/.netrc and is never stored here.

if (file.exists("DESCRIPTION") && requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(rGEDI)
}

earthdata_login()

granule_ids <- sprintf(
  "GEDI02_A_2019108080339_O01964_%02d_T05337_02_003_01_V002",
  1:4
)
urls <- sprintf(
  paste0(
    "https://data.lpdaac.earthdatacloud.nasa.gov/lp-prod-protected/",
    "GEDI02_A.002/%s/%s.h5"
  ),
  granule_ids, granule_ids
)

cache_dir <- file.path("readme", ".orbit-cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

parts <- lapply(seq_along(urls), function(i) {
  cache_file <- file.path(cache_dir, sprintf("part-%02d.csv", i))
  if (file.exists(cache_file)) {
    message("Using cached coordinate sample for orbit part ", i, " of 4")
    return(data.table::fread(cache_file))
  }
  for (attempt in 1:4) {
    message(
      "Streaming orbit part ", i, " of ", length(urls), ": ",
      granule_ids[i], " (attempt ", attempt, ")"
    )
    h5 <- NULL
    track <- tryCatch({
      h5 <- readLevel2A(urls[i])
      getGEDITrack(h5, every = 200)
    }, error = function(e) {
      message("Earthdata request failed: ", conditionMessage(e))
      NULL
    }, finally = {
      if (!is.null(h5)) try(close(h5), silent = TRUE)
    })
    if (!is.null(track)) {
      track$granule_part <- i
      track$granule_id <- granule_ids[i]
      track$track <- track$beam
      data.table::fwrite(track, cache_file)
      return(track)
    }
    if (attempt < 4) Sys.sleep(5 * attempt)
  }
  stop("Could not stream orbit part ", i, " after four attempts.")
})

complete_orbit <- data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
complete_orbit[, track := beam]
data.table::setorder(complete_orbit, track, sequence)
data.table::fwrite(complete_orbit, "readme/gedi-orbit-track.csv")
unlink(cache_dir, recursive = TRUE)

plot_gedi_orbit_animation(
  complete_orbit,
  output_file = "readme/gedi-orbit-animation.gif",
  title = "GEDI aboard the International Space Station",
  duration = 12,
  earth_rotation_speed = 2,
  launch = FALSE
)
plot_gedi_orbit_animation(
  complete_orbit,
  output_file = "readme/gedi-orbit-animation.html",
  title = "GEDI aboard the International Space Station",
  duration = 12,
  track_speed = 2,
  earth_rotation_speed = 2,
  launch = FALSE
)

message(
  "Complete orbit built from ", nrow(complete_orbit), " sampled footprints; ",
  "latitude ", round(min(complete_orbit$latitude), 3), " to ",
  round(max(complete_orbit$latitude), 3), " degrees."
)

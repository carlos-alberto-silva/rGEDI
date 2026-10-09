# CRAN-safe GEDI waveform simulator -------------------------------------

.read_als_points <- function(input) {
  if (inherits(input, "LAS")) return(as.data.frame(input@data))
  if (inherits(input, c("data.frame", "data.table"))) return(as.data.frame(input))
  if (inherits(input, "SpatVector")) return(cbind(terra::crds(input), terra::values(input)))
  if (is.character(input)) {
    if (!all(file.exists(input))) stop("ALS input file does not exist.")
    if (!requireNamespace("lidR", quietly = TRUE)) stop("Install 'lidR' to read LAS/LAZ files.")
    clouds <- lapply(input, lidR::readLAS)
    return(do.call(rbind, lapply(clouds, function(x) as.data.frame(x@data))))
  }
  stop("'input' must be LAS/LAZ paths, a lidR LAS object, SpatVector, or table.")
}

.als_xyz <- function(points) {
  nms <- tolower(names(points))
  find <- function(options) {
    hit <- match(options, nms, nomatch = 0L)
    hit <- hit[hit > 0L]
    if (length(hit)) names(points)[hit[[1L]]] else NA_character_
  }
  cols <- c(find(c("x", "easting", "lon", "longitude")),
            find(c("y", "northing", "lat", "latitude")),
            find(c("z", "elevation", "height")))
  if (anyNA(cols)) stop("ALS data must contain X, Y, and Z columns.")
  cols
}

.waveform_centres <- function(points, xyz, coords, listCoord, gridBound, gridStep) {
  if (!is.null(coords)) return(matrix(as.numeric(coords), ncol = 2L, byrow = TRUE))
  if (!is.null(listCoord)) return(as.matrix(utils::read.table(listCoord)[, 1:2, drop = FALSE]))
  if (!is.null(gridBound)) {
    xs <- seq(gridBound[1L], gridBound[2L], by = gridStep)
    ys <- seq(gridBound[3L], gridBound[4L], by = gridStep)
    return(as.matrix(expand.grid(x = xs, y = ys)))
  }
  matrix(c(mean(points[[xyz[1L]]], na.rm = TRUE), mean(points[[xyz[2L]]], na.rm = TRUE)), ncol = 2L)
}

.gaussian_kernel <- function(sigma_bins) {
  if (!is.finite(sigma_bins) || sigma_bins <= 0) return(1)
  radius <- max(1L, ceiling(4 * sigma_bins))
  k <- stats::dnorm(seq.int(-radius, radius), sd = sigma_bins)
  k / sum(k)
}

.simulate_one_waveform <- function(points, xyz, centre, fSigma, res, pulse_sigma,
                                   decimate, maxBins, countOnly) {
  dx <- points[[xyz[1L]]] - centre[1L]
  dy <- points[[xyz[2L]]] - centre[2L]
  keep <- is.finite(dx) & is.finite(dy) & is.finite(points[[xyz[3L]]]) &
    (dx * dx + dy * dy <= (4 * fSigma)^2)
  p <- points[keep, , drop = FALSE]
  if (!nrow(p)) return(NULL)
  if (decimate < 1) p <- p[stats::runif(nrow(p)) <= decimate, , drop = FALSE]
  if (!nrow(p)) return(NULL)
  dx <- p[[xyz[1L]]] - centre[1L]
  dy <- p[[xyz[2L]]] - centre[2L]
  weight <- if (countOnly) rep(1, nrow(p)) else exp(-0.5 * (dx * dx + dy * dy) / fSigma^2)
  z <- p[[xyz[3L]]]
  low <- floor(min(z) / res) * res
  high <- ceiling(max(z) / res) * res
  breaks <- seq(low, high + res, by = res)
  if (length(breaks) - 1L > maxBins) breaks <- seq(low, high, length.out = maxBins + 1L)
  bin <- pmax(1L, pmin(length(breaks) - 1L, findInterval(z, breaks, all.inside = TRUE)))
  amplitude <- numeric(length(breaks) - 1L)
  for (i in seq_along(bin)) amplitude[bin[i]] <- amplitude[bin[i]] + weight[i]
  amplitude <- stats::filter(amplitude, .gaussian_kernel(pulse_sigma / res), sides = 2, circular = FALSE)
  amplitude[is.na(amplitude)] <- 0
  elevation <- breaks[-length(breaks)] + diff(breaks) / 2
  list(amplitude = as.numeric(amplitude), elevation = elevation,
       elevation_bin0 = max(elevation), elevation_lastbin = min(elevation))
}

.write_simulated_l1b <- function(waves, centres, output) {
  h5 <- hdf5r::H5File$new(output, mode = "w")
  beam <- h5$create_group("BEAM0000")
  geo <- beam$create_group("geolocation")
  counts <- as.integer(vapply(waves, function(x) length(x$amplitude), integer(1)))
  starts <- cumsum(c(1L, head(counts, -1L)))
  beam[["shot_number"]] <- seq.int(0, length(waves) - 1L)
  beam[["rx_sample_count"]] <- counts
  beam[["rx_sample_start_index"]] <- starts
  beam[["rxwaveform"]] <- unlist(lapply(waves, `[[`, "amplitude"), use.names = FALSE)
  geo[["elevation_bin0"]] <- vapply(waves, `[[`, numeric(1), "elevation_bin0")
  geo[["elevation_lastbin"]] <- vapply(waves, `[[`, numeric(1), "elevation_lastbin")
  geo[["longitude_bin0"]] <- centres[, 1L]
  geo[["latitude_bin0"]] <- centres[, 2L]
  h5$close_all()
  new("gedi.level1b", h5 = hdf5r::H5File$new(output, mode = "r"))
}

#' Simulate GEDI-like full waveforms from airborne laser scanning data
#'
#' This portable implementation uses Gaussian footprint weighting and pulse
#' convolution. It restores waveform simulation without bundling the original
#' native simulator's GSL, GDAL, GeoTIFF, HDF5, SZIP, and zlib build stack.
#'
#' @param input LAS/LAZ paths, a `lidR::LAS`, `SpatVector`, or table with X/Y/Z.
#' @param output Output HDF5 path.
#' @param ground Retained for compatibility; ground returns remain in the waveform.
#' @param ascii Return/write a simple ASCII waveform instead of HDF5.
#' @param waveID Optional waveform identifiers.
#' @param coords A coordinate pair or matrix of footprint centers.
#' @param listCoord Text file containing footprint centers.
#' @param gridBound Optional `c(xmin, xmax, ymin, ymax)` simulation grid.
#' @param gridStep Grid spacing in input coordinate units.
#' @param pSigma Pulse sigma in meters. A negative value derives it from `pFWHM`.
#' @param pFWHM Pulse full width at half maximum in nanoseconds.
#' @param fSigma Footprint sigma in input coordinate units.
#' @param res Vertical waveform resolution.
#' @param decimate Probability of retaining an ALS return.
#' @param maxBins Maximum number of waveform bins.
#' @param countOnly Use unweighted return counts inside the footprint.
#' @param keepOld Do not overwrite an existing output file.
#' @param seed Optional random seed.
#' @param ... Compatibility arguments accepted from the original simulator.
#' @return A [`gedi.level1b-class`] object, or a waveform table when `ascii=TRUE`.
#' @references Hancock et al. (2019) \doi{10.1029/2018EA000506}
#' @export
gediWFSimulator <- function(input, output = tempfile(fileext = ".h5"), ground = TRUE,
                            ascii = FALSE, waveID = NULL, coords = NULL,
                            listCoord = NULL, gridBound = NULL, gridStep = 30,
                            pSigma = -1, pFWHM = 15, fSigma = 5.5, res = 0.15,
                            decimate = 1, maxBins = 1024L, countOnly = FALSE,
                            keepOld = FALSE, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  if (!isTRUE(ascii) && !grepl("\\.h5$", output, ignore.case = TRUE)) output <- paste0(output, ".h5")
  if (file.exists(output) && keepOld) return(.read_local_gedi(output, "GEDI01_B"))
  if (file.exists(output)) unlink(output)
  if (!is.finite(fSigma) || fSigma <= 0 || !is.finite(res) || res <= 0) {
    stop("`fSigma` and `res` must be positive finite numbers.")
  }
  if (!is.finite(decimate) || decimate <= 0 || decimate > 1) {
    stop("`decimate` must be greater than zero and no greater than one.")
  }
  points <- .read_als_points(input)
  xyz <- .als_xyz(points)
  centres <- .waveform_centres(points, xyz, coords, listCoord, gridBound, gridStep)
  if (ncol(centres) != 2L) stop("Footprint coordinates must have two columns.")
  pulse_sigma <- if (pSigma > 0) pSigma else pFWHM * 0.15 / 2.355
  waves <- lapply(seq_len(nrow(centres)), function(i) {
    .simulate_one_waveform(points, xyz, centres[i, ], fSigma, res, pulse_sigma, decimate, maxBins, countOnly)
  })
  valid <- !vapply(waves, is.null, logical(1))
  waves <- waves[valid]
  centres <- centres[valid, , drop = FALSE]
  if (!length(waves)) stop("No ALS returns intersect the requested GEDI footprints.")
  if (isTRUE(ascii)) {
    ans <- data.table::rbindlist(lapply(seq_along(waves), function(i) data.table::data.table(
      waveID = if (is.null(waveID)) i - 1L else waveID[i],
      x = centres[i, 1L], y = centres[i, 2L], elevation = waves[[i]]$elevation,
      amplitude = waves[[i]]$amplitude)))
    utils::write.table(ans, output, row.names = FALSE)
    return(ans)
  }
  .write_simulated_l1b(waves, centres, output)
}

.weighted_quantile <- function(x, w, probs) {
  o <- order(x)
  x <- x[o]; w <- pmax(0, w[o])
  if (!sum(w)) return(rep(NA_real_, length(probs)))
  stats::approx(cumsum(w) / sum(w), x, xout = probs, rule = 2, ties = "ordered")$y
}

#' Compute metrics from GEDI or simulated full waveforms
#'
#' @param input A [`gedi.level1b-class`] object or list of such objects.
#' @param outRoot Optional output filename root.
#' @param rhRes Relative-height percentile spacing.
#' @param ground Logical; include the estimated ground elevation.
#' @param rhoG,rhoC Ground and canopy reflectance.
#' @param ... Compatibility arguments from the original native implementation.
#' @return A [data.table::data.table] with waveform metrics.
#' @export
gediWFMetrics <- function(input, outRoot = NULL, rhRes = 5, ground = TRUE,
                          rhoG = 0.4, rhoC = 0.57, ...) {
  objects <- if (is(input, "gedi.level1b")) list(input) else input
  if (!is.list(objects) || !all(vapply(objects, is, logical(1), "gedi.level1b"))) stop("'input' must contain GEDI Level 1B objects.")
  probs <- seq(0, 100, by = rhRes)
  rows <- list(); k <- 0L
  for (obj in objects) {
    h5 <- obj@h5
    for (beam in .gedi_beams(h5)) {
      shots <- h5[[paste0(beam, "/shot_number")]][]
      for (shot in shots) {
        wf <- getLevel1BWF(obj, shot)@dt
        amp <- pmax(0, wf$rxwaveform)
        elev <- wf$elevation
        q <- .weighted_quantile(elev, amp, probs / 100)
        ground_elev <- q[[1L]]
        canopy <- elev > ground_elev + 2
        e_can <- sum(amp[canopy]); e_ground <- sum(amp[!canopy])
        cover <- e_can / (e_can + e_ground * rhoC / rhoG)
        p <- amp / sum(amp)
        fhd <- -sum(p[p > 0] * log(p[p > 0]))
        k <- k + 1L
        row <- data.table::data.table(beam = beam, shot_number = shot,
          gHeight = if (ground) ground_elev else NA_real_,
          signal_top = max(elev[amp > 0]), signal_bottom = min(elev[amp > 0]),
          cover = cover, waveEnergy = sum(amp), FHD = fhd)
        for (i in seq_along(probs)) row[[paste0("rh", probs[i])]] <- q[i] - ground_elev
        rows[[k]] <- row
      }
    }
  }
  ans <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (!is.null(outRoot)) utils::write.table(ans, paste0(outRoot, ".metric.txt"), row.names = FALSE)
  ans
}

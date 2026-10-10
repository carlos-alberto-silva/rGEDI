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

.als_field <- function(points, options) {
  nms <- tolower(names(points))
  hit <- match(tolower(options), nms, nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit)) names(points)[hit[[1L]]] else NULL
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

.weighted_bins <- function(bin, weight, nbins) {
  out <- numeric(nbins)
  if (!length(bin)) return(out)
  values <- rowsum(weight, bin, reorder = FALSE)
  out[as.integer(rownames(values))] <- values[, 1L]
  out
}

.simulate_one_waveform <- function(points, xyz, centre, fSigma, res, pulse_sigma,
                                   decimate, maxBins, method, normalize,
                                   density_res, intensity_threshold, ground) {
  dx <- points[[xyz[1L]]] - centre[1L]
  dy <- points[[xyz[2L]]] - centre[2L]
  keep <- is.finite(dx) & is.finite(dy) & is.finite(points[[xyz[3L]]]) &
    exp(-0.5 * (dx * dx + dy * dy) / fSigma^2) >= intensity_threshold
  p <- points[keep, , drop = FALSE]
  if (!nrow(p)) return(NULL)
  if (decimate < 1) p <- p[stats::runif(nrow(p)) <= decimate, , drop = FALSE]
  if (!nrow(p)) return(NULL)
  dx <- p[[xyz[1L]]] - centre[1L]
  dy <- p[[xyz[2L]]] - centre[2L]
  footprint <- exp(-0.5 * (dx * dx + dy * dy) / fSigma^2)
  weight <- footprint
  if (method == "intensity") {
    intensity <- .als_field(p, c("intensity", "refl", "reflectance"))
    if (is.null(intensity)) stop("Intensity simulation requires an ALS intensity column.")
    weight <- weight * as.numeric(p[[intensity]])
  } else if (method == "fraction") {
    nreturns <- .als_field(p, c("numberofreturns", "number_of_returns", "nret"))
    if (!is.null(nreturns)) weight <- weight / pmax(1, as.numeric(p[[nreturns]]))
  }
  if (isTRUE(normalize)) {
    gx <- floor((p[[xyz[1L]]] - min(p[[xyz[1L]]])) / density_res)
    gy <- floor((p[[xyz[2L]]] - min(p[[xyz[2L]]])) / density_res)
    cell <- interaction(gx, gy, drop = TRUE)
    return_number <- .als_field(p, c("returnnumber", "return_number", "retnumb"))
    number_returns <- .als_field(p, c("numberofreturns", "number_of_returns", "nret"))
    beam_end <- if (!is.null(return_number) && !is.null(number_returns)) {
      as.integer(p[[return_number]]) == as.integer(p[[number_returns]])
    } else rep(TRUE, nrow(p))
    density <- tabulate(as.integer(cell)[beam_end], nbins = nlevels(cell))
    density[density < 1L] <- 1L
    weight <- weight / density[as.integer(cell)]
  }
  z <- p[[xyz[3L]]]
  low <- floor((min(z) - 35) / res) * res
  high <- ceiling((max(z) + 35) / res) * res
  breaks <- seq(low, high + res, by = res)
  bin <- pmax(1L, pmin(length(breaks) - 1L, findInterval(z, breaks, all.inside = TRUE)))
  amplitude <- .weighted_bins(bin, weight, length(breaks) - 1L)
  class_col <- .als_field(p, c("classification", "class"))
  ground_point <- if (!isTRUE(ground) || is.null(class_col)) {
    rep(FALSE, nrow(p))
  } else as.integer(p[[class_col]]) == 2L
  ground_amplitude <- .weighted_bins(
    bin[ground_point], weight[ground_point], length(breaks) - 1L
  )
  amplitude <- stats::filter(amplitude, .gaussian_kernel(pulse_sigma / res), sides = 2, circular = FALSE)
  ground_amplitude <- stats::filter(ground_amplitude, .gaussian_kernel(pulse_sigma / res), sides = 2, circular = FALSE)
  amplitude[is.na(amplitude)] <- 0
  ground_amplitude[is.na(ground_amplitude)] <- 0
  elevation <- breaks[-length(breaks)] + res / 2
  if (length(amplitude) > maxBins) {
    active <- which(amplitude > max(amplitude) * 1e-10)
    centre_bin <- if (length(active)) {
      as.integer(round(mean(range(active))))
    } else {
      as.integer(ceiling(length(amplitude) / 2))
    }
    first <- max(1L, min(length(amplitude) - maxBins + 1L,
                         centre_bin - as.integer(floor(maxBins / 2))))
    take <- seq.int(first, length.out = maxBins)
    amplitude <- amplitude[take]
    ground_amplitude <- ground_amplitude[take]
    elevation <- elevation[take]
  }
  integral <- sum(amplitude) * res
  if (integral > 0) {
    amplitude <- amplitude / integral
    ground_amplitude <- ground_amplitude / integral
  }
  # GEDI L1B stores the upper RX-window edge in elevation_bin0; the first
  # sample centre is one interval below it and samples proceed downward.
  list(amplitude = rev(as.numeric(amplitude)), ground = rev(as.numeric(ground_amplitude)),
       elevation = rev(elevation),
       elevation_bin0 = max(elevation) + res,
       elevation_lastbin = min(elevation))
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
  beam[["grxwaveform"]] <- unlist(lapply(waves, `[[`, "ground"), use.names = FALSE)
  geo[["elevation_bin0"]] <- vapply(waves, `[[`, numeric(1), "elevation_bin0")
  geo[["elevation_lastbin"]] <- vapply(waves, `[[`, numeric(1), "elevation_lastbin")
  geo[["longitude_bin0"]] <- centres[, 1L]
  geo[["latitude_bin0"]] <- centres[, 2L]
  h5[["FSIGMA"]] <- attr(waves, "fSigma")
  h5[["PSIGMA"]] <- attr(waves, "pSigma")
  h5[["WAVEFORM_RES"]] <- attr(waves, "res")
  h5$close_all()
  new("gedi.level1b", h5 = hdf5r::H5File$new(output, mode = "r"))
}

#' Simulate GEDI-like full waveforms from airborne laser scanning data
#'
#' This portable implementation follows the core `gediRat` algorithm: Gaussian
#' footprint weighting, optional return-density normalization, vertical
#' binning, pulse convolution, separation of LAS class 2 ground returns, and
#' unit-integral normalization. Vectorized binning replaces the original point
#' loop while avoiding its native GSL/GDAL build stack.
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
#' @param countOnly Use count-based returns rather than ALS intensity. Gaussian
#'   footprint weights are still applied.
#' @param method Waveform contribution: `"count"`, `"intensity"`, or
#'   `"fraction"`. `countOnly = TRUE` selects `"count"` for compatibility.
#' @param normalize Correct for local ALS sampling density, as in the original
#'   simulator's default behavior.
#' @param density_res Cell size used for sampling-density normalization.
#' @param intensity_threshold Smallest retained Gaussian footprint weight.
#' @param noise Relative Gaussian noise standard deviation. For example,
#'   `noise = 0.03` adds noise with a standard deviation equal to three percent
#'   of the waveform maximum. Zero leaves the waveform unchanged.
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
                            method = c("count", "intensity", "fraction"),
                            normalize = TRUE, density_res = 1,
                            intensity_threshold = 0.0006,
                            noise = 0, keepOld = FALSE, seed = NULL, ...) {
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
  if (!is.finite(density_res) || density_res <= 0) {
    stop("`density_res` must be a positive finite number.")
  }
  if (!is.finite(intensity_threshold) || intensity_threshold <= 0 ||
      intensity_threshold > 1) {
    stop("`intensity_threshold` must be greater than zero and no greater than one.")
  }
  maxBins <- as.integer(maxBins)[1L]
  if (is.na(maxBins) || maxBins < 2L) stop("`maxBins` must be an integer of at least two.")
  noise <- as.numeric(noise)[1L]
  if (!is.finite(noise) || noise < 0) stop("`noise` must be a non-negative number.")
  points <- .read_als_points(input)
  method <- match.arg(method)
  if (isTRUE(countOnly)) method <- "count"
  xyz <- .als_xyz(points)
  centres <- .waveform_centres(points, xyz, coords, listCoord, gridBound, gridStep)
  if (ncol(centres) != 2L) stop("Footprint coordinates must have two columns.")
  pulse_sigma <- if (pSigma > 0) pSigma else pFWHM * 0.15 / 2.355
  if (!is.finite(pulse_sigma) || pulse_sigma <= 0) {
    stop("`pSigma`, or the value derived from `pFWHM`, must be positive.")
  }
  waves <- lapply(seq_len(nrow(centres)), function(i) {
    .simulate_one_waveform(
      points, xyz, centres[i, ], fSigma, res, pulse_sigma, decimate,
      maxBins, method, normalize, density_res, intensity_threshold, ground
    )
  })
  valid <- !vapply(waves, is.null, logical(1))
  waves <- waves[valid]
  centres <- centres[valid, , drop = FALSE]
  if (!length(waves)) stop("No ALS returns intersect the requested GEDI footprints.")
  if (noise > 0) {
    waves <- lapply(waves, function(wave) {
      sigma <- noise * max(wave$amplitude, na.rm = TRUE)
      wave$amplitude <- wave$amplitude + stats::rnorm(length(wave$amplitude), sd = sigma)
      wave
    })
  }
  attr(waves, "fSigma") <- fSigma
  attr(waves, "pSigma") <- pulse_sigma
  attr(waves, "res") <- res
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
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]; w <- w[keep]
  o <- order(x)
  x <- x[o]; w <- pmax(0, w[o])
  if (!sum(w)) return(rep(NA_real_, length(probs)))
  stats::approx(cumsum(w) / sum(w), x, xout = probs, rule = 2, ties = "ordered")$y
}

.wave_ground <- function(h5, beam, index, elevation, amplitude, method) {
  count <- as.integer(h5[[paste0(beam, "/rx_sample_count")]][index])
  start <- as.integer(h5[[paste0(beam, "/rx_sample_start_index")]][index])
  path <- paste0(beam, "/grxwaveform")
  truth <- NULL
  if (method %in% c("auto", "truth") && h5$exists(path)) {
    truth <- h5[[path]][start:(start + count - 1L)]
    if (length(truth) && sum(pmax(0, truth)) > 0) {
      return(list(height = stats::weighted.mean(elevation, pmax(0, truth)),
                  ground = truth, source = "classified_ground"))
    }
  }
  positive <- pmax(0, amplitude)
  if (method == "quantile") {
    return(list(height = .weighted_quantile(elevation, positive, 0.01),
                ground = truth, source = "quantile"))
  }
  smooth <- stats::filter(positive, .gaussian_kernel(max(1, 0.57 / abs(diff(elevation)[1]))), sides = 2)
  smooth[is.na(smooth)] <- 0
  peaks <- which(smooth >= c(-Inf, head(smooth, -1L)) &
                   smooth > c(tail(smooth, -1L), -Inf) &
                   smooth >= max(smooth) * 0.005)
  height <- if (length(peaks)) min(elevation[peaks]) else .weighted_quantile(elevation, positive, 0.01)
  list(height = height, ground = truth, source = "lowest_peak")
}

#' Compute metrics from GEDI or simulated full waveforms
#'
#' @param input A [`gedi.level1b-class`] object or list of such objects.
#' @param outRoot Optional output filename root.
#' @param rhRes Relative-height percentile spacing.
#' @param ground Logical; include the estimated ground elevation.
#' @param rhoG,rhoC Ground and canopy reflectance.
#' @param ground_method Ground estimator. `"auto"` uses the classified ground
#'   waveform when present and otherwise the lowest smoothed waveform peak.
#' @param ... Compatibility arguments from the original native implementation.
#' @return A [data.table::data.table] with waveform metrics.
#' @export
gediWFMetrics <- function(input, outRoot = NULL, rhRes = 5, ground = TRUE,
                          rhoG = 0.4, rhoC = 0.57,
                          ground_method = c("auto", "truth", "lowest_peak", "quantile"), ...) {
  objects <- if (is(input, "gedi.level1b")) list(input) else input
  if (!is.list(objects) || !all(vapply(objects, is, logical(1), "gedi.level1b"))) stop("'input' must contain GEDI Level 1B objects.")
  probs <- seq(0, 100, by = rhRes)
  ground_method <- match.arg(ground_method)
  rows <- list(); k <- 0L
  for (obj in objects) {
    h5 <- obj@h5
    for (beam in .gedi_beams(h5)) {
      shots <- h5[[paste0(beam, "/shot_number")]][]
      for (shot_index in seq_along(shots)) {
        shot <- shots[[shot_index]]
        wf <- getLevel1BWF(obj, shot)@dt
        amp <- pmax(0, wf$rxwaveform)
        elev <- wf$elevation
        ground_info <- .wave_ground(h5, beam, shot_index, elev, amp, ground_method)
        if (ground_method == "truth" && ground_info$source != "classified_ground") {
          stop("No classified-ground waveform is available for shot ", shot, ".")
        }
        ground_elev <- ground_info$height
        q <- .weighted_quantile(elev, amp, probs / 100)
        rh <- q - ground_elev
        rh[rh < 0] <- 0
        if (!is.null(ground_info$ground) && sum(pmax(0, ground_info$ground)) > 0) {
          e_ground <- sum(pmax(0, ground_info$ground))
          e_can <- sum(amp) - e_ground
        } else {
          ground_bin <- which.min(abs(elev - ground_elev))
          e_ground <- min(sum(amp), 2 * sum(amp[seq_len(ground_bin)]))
          e_can <- sum(amp) - e_ground
        }
        e_can <- max(0, e_can); e_ground <- max(0, e_ground)
        cover <- e_can / (e_can + e_ground * rhoC / rhoG)
        p <- amp / sum(amp)
        fhd <- -sum(p[p > 0] * log(p[p > 0]))
        k <- k + 1L
        row <- data.table::data.table(beam = beam, shot_number = shot,
          gHeight = if (ground) ground_elev else NA_real_,
          ground_method = ground_info$source,
          signal_top = max(elev[amp > 0]), signal_bottom = min(elev[amp > 0]),
          cover = cover, waveEnergy = sum(amp) * abs(diff(elev)[1]), FHD = fhd)
        for (i in seq_along(probs)) row[[paste0("rh", probs[i])]] <- rh[i]
        rows[[k]] <- row
      }
    }
  }
  ans <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (!is.null(outRoot)) utils::write.table(ans, paste0(outRoot, ".metric.txt"), row.names = FALSE)
  ans
}

#'GEDI full-waveform data simulation
#'
#'@description Simulate GEDI full-waveform data from Airborne Laser Scanning (ALS) 3-D point cloud
#'
#' Input output filenames and format
#' @param input character vector. lasfile input filename
# inList should be parsed from vector input
# @param inList list. input file list (ASCII file) for multiple files
#' @param output character. output filename

#Ground should always be true
# @param ground record separate ground and canopy waveforms

#Always HDF
# @param hdf write output as HDF5. Best with gridded or list of coords
# @param ascii write output as ASCII (default). Good for quick tests
#' @param waveID id. supply a waveID to pass to the output (only for single footprints)
#' Single footprint, list of footprints, or grid of footprints
#' @param coord lon lat numeric vector. footprint coordinate in same system as lasfile
#' @param listCoord name. Text file with list of coordinates. Pattern: X Y `[waveID]` `[geoCoordsX]` `[geoCoordsY]`. `[]` are optional, separated by spaces.
#' @param gridBound minX maxX minY maxY numeric vector. make a grid of waveforms in this box
#' @param gridStep res. grid step size
#' Lidar characteristics. Defaults are expected GEDI values. pSigmasig. set Gaussian pulse width as 1 sigma
#' @param pFWHM fhwm. set Gaussian pulse width as FWHM in ns
#' @param readPulse file. read pulse shape and width from a file insteda of making Gaussian
#' @param fSigma sig. set footprint width
#' @param wavefront file. read wavefront shape from file instead of setting Gaussian. Note that footprint width is still set by fSigma
#' @param res res. range resolution of waveform digitisation to output, in units of ALS data

# Not LVIS
# @param LVIS use LVIS pulse length, sigma=6. 25m
#' @param topHat use a top hat wavefront
#' @param sideLobe use side lobes
#' @param lobeAng ang. lobe axis azimuth
#' Input data quality filters
#' @param checkCover check that the footprint is covered by ALS data. Do not output if not
#' @param maxScanAng ang. maximum scan angle, degrees
#' @param decimate x. probability of accepting an ALS beam
#' Computational speed options
#' @param pBuff s. point reading buffer size in Gbytes
#' @param maxBins for HDF5, limit number of bins to save trimming.
#' @param countOnly only use count method
#' @param pulseAfter apply the pulse smoothing after binning for computational speed, at the risk of aliasing (default)
#' @param pulseBefore apply the pulse smoothing before binning to avoid the risk of aliasing, at the expense of computational speed
#' @param noNorm don't normalise for ALS density

#' Octree
#' @param noOctree do not use an octree
#' @param octLevels n. number of octree levels to use
#' @param nOctPix n. number of octree pixels along a side for the top level

#' Using full-waveform input data (not tested)
# Not supported yet
# @param decon deconvolve
# @param indDecon deconvolve individual beams
# @param readWave read full-waveform where available
# Miscellaneous

# R user should never only list files
# @param listFiles list files. Do not read them
#' @param keepOld do not overwrite old files, if they exist
#' @param useShadow account for shadowing in discrete return data through voxelisation
#' @param polyGround find mean ground elevation and slope through fitting a polynomial

# nnGround is not working yet
# @param nnGround find mean ground elevation and slope through nearest neighbour
#' @param seed n integer. random number seed
#'
#' #'
#' @return A S4 object of class \code{\link[hdf5r]{hdf5rfile}} in the \emph{hdf5r} package.
#'
#' @seealso
#' i) Hancock, S., Armston, J., Hofton, M., Sun, X., Tang, H., Duncanson, L.I., Kellner,
#' J.R. and Dubayah, R., 2019. The GEDI simulator: A large‐footprint waveform lidar simulator
#' for calibration and validation of spaceborne missions. Earth and Space Science.
#' https://doi.org/10.1029/2018EA000506
#'
#' ii) gediSimulator: https://bitbucket.org/StevenHancock/gedisimulator/src/master/
#'
#' @examples
#'# specify the path to ALS data
#'LASfile <- system.file("extdata", "LASexample1.las", package="rGEDI")
#'
#'Reading and plot LASfile
#'library(lidR)
#'LAS<-readLAS(LASfile)
#'plot(LAS)
#'
#'# Simulate GEDI full-waveform
#'wf<-gediWFSimulator(input=LASfile,output="gediSimulation.h5")
#'
#' @import hdf5r
#' @useDynLib rGEDI
#' @export
gediWFSimulator = function(
  input,
  output,
  waveID = NULL,
  coords = NULL,
  listCoord = NULL,
  gridBound = NULL,
  gridStep = 30.0,

  pSigma = -1.0,
  pFWHM = 15.0,
  readPulse = NULL,
  fSigma = 5.5,
  wavefront = NULL,
  res = 0.15,
  topHat = FALSE,
  sideLobe = FALSE,
  lobeAng = 0.0,

  checkCover = FALSE,
  maxScanAng = 1000000.0,
  decimate = 1.0,

  pBuff = as.integer(200000000),
  maxBins = as.integer(1024),
  countOnly = FALSE,
  pulseAfter = FALSE,
  pulseBefore = TRUE,
  noNorm = FALSE,

  noOctree = FALSE,
  octLevels = as.integer(0),
  nOctPix = as.integer(40),

  keepOld = FALSE,
  useShadow = FALSE,
  polyGround = FALSE,
  seed = NULL) {

  # Set parameters that shouldn't be changed for GEDI
  ground = TRUE
  hdf = TRUE
  ascii = FALSE
  LVIS = FALSE

  decon = FALSE
  indDecon = FALSE
  readWave = FALSE

  listFiles = FALSE
  nnGround = FALSE

  # Check values
  stopifnotMessage(
    all(file.exists(input)),
    all(fs::path_ext(input) == "las"),
    dir.exists(fs::path_dir(output)),
    is.null(waveID) || length(coords) == 2, # If waveID should only work along with coords
    checkNumericLength(coords, 2),
    is.null(listCoord) || file.exists(listCoord),
    checkNumericLength(gridBound, 4),
    checkNumeric(gridStep),
    checkNumeric(pFWHM),
    checkFilepath(readPulse, newFile=FALSE, optional=TRUE),
    checkNumeric(fSigma),
    checkFilepath(wavefront, newFile=FALSE, optional=TRUE),
    checkNumeric(res),
    checkLogical(topHat),
    checkLogical(sideLobe),
    checkNumeric(lobeAng),
    checkLogical(checkCover),
    checkNumeric(maxScanAng),
    checkNumeric(decimate),
    checkNumeric(pBuff),
    checkInteger(maxBins),
    checkLogical(countOnly),
    checkLogical(pulseAfter),
    checkLogical(pulseBefore),
    checkLogical(noNorm),
    checkLogical(noOctree),
    checkInteger(octLevels),
    checkInteger(nOctPix),
    checkLogical(keepOld),
    checkLogical(useShadow),
    checkLogical(polyGround),
    checkInteger(seed)
  )

  if (is.null(coords) && is.null(listCoord) && is.null(gridBound)) {
    stop("Coordinates for the waveforms should be provided!\nTIP: Use coords, listCoord or gridBound.")
  }

  inputInList = inputOrInList(input)
  if (fs::path_ext(output) != "h5") {
    output = paste0(output, ".h5")
  }

  .Call("C_gediSimulator",
        inputInList[[1]],
        output,
        inputInList[[2]],
        ground,
        hdf,
        ascii,
        waveID,

        coords,
        listCoord,
        gridBound,
        gridStep,

        pSigma,
        pFWHM,
        readPulse,
        fSigma,
        wavefront,
        res,
        LVIS,
        topHat,
        sideLobe,
        lobeAng,

        checkCover,
        maxScanAng,
        decimate,

        as.integer(pBuff),
        as.integer(maxBins),
        countOnly,
        pulseAfter,
        pulseBefore,
        noNorm,

        noOctree,
        as.integer(octLevels),
        as.integer(nOctPix),

        decon,
        indDecon,
        readWave,

        listFiles,
        keepOld,
        useShadow,
        polyGround,
        nnGround,
        as.integer(seed))

  unloadLibrary()

  cleanInList(inputInList)
  result = tryCatch(hdf5r::H5File$new(output, "r+"), error=function(e) stop("The output file was not created\nSomething went wrong!"))

  result<- new("gedi.level1bSim", h5 = result)

  return(result)
}

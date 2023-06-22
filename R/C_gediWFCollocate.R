#' Collocate waveforms based in ALS
#'
#' @description Uses the correlation method in Blair and Hofton (1999) to colocate a large-footprint
#' lidar dataset with a small-footprint, discrete-return dataset. Note that it requires the full-waveform
#' LVIS or GEDI data, which is contained in the L1B files. It uses the Pearson correlation to find the best
#' affine transformation (x and y only, or x, y and z) and footprint size needed to align a large-footprint
#' dataset with a small-footprint dataset. It has three potential modes of operation.
#'
#'  - It can test a grid of affine transformations and give the correlation for every point
#'(for a single footprint width), as used in Blair and Hofton (1999).
#'  - It can use a simplex to move along the error surface and find the optimum transformation
#' and footprint size. Note that initial location needs to be within around 20 m of the true location
#' for their to be a sufficient gradient on the error surface.
#'  - It can use a hybrid of the two, testing every step on a coarse grid, then setting a simplex off
#' from the location of maximum correlation. This allows a rapid assessment to get within 20 m of the
#' true value, then uses the simplex to find the offsets precisely.
#'  - It requires the large-footprint system pulse shape (either a Gaussian width or a file containing
#' range and intensity) and the EPSG codes for the two datasets. It reads ALS data in .las format and
#' can read either simulated HDF5 files from gediRat or LVIS in HDF5 or lgw format. A reader for GEDI
#' data will be added once that data is available.
#'
#' If the full grid is used, it outputs an ASCII file with the correlation for each x, y and z offset.
#' If it uses a simplex it outputs the single optimum offset.
#'
#  Input-output
#' @param output           output filename
#' @param als file;        input als file
#' @param gedi file;       single input GEDI/LVIS L1B file
#' @param readHDFgedi;     read GEDI HDF5 input (default)
#' @param lgw;             LVIS is in lgw (default is GEDI hdf5)
#' @param readHDFlvis;     read GEDI HDF5 input (default is GEDI hdf5)
#' @param bounds minX minY maxX maxY;    bounds to use, specified in ALS projection
#' @param leaveEmpty;      exit if there are no usable footprints
#' @param lEPSG epsg;      LVIS projection
#' @param aEPSG epsg;      ALS projection
#' @param beamList 11111111; 0/1 for whether or not to use beams 1-8
#' @param skipBeams n;     list of beam numbers to skip. No spaces between (eg 123)
#' @param readBeams n;     list of beam numbers to read. No spaces between (eg 123)
#'
#  Multi-mode options
#' @param solveCofG        use centre of gravity to match vertical offsets
#' @param maxVshift x;     grid or geoError mode, vertical distance to search over
#' @param vStep z;         grid or geoError mode, vertical step size
#'
#  Grid mode operation (multi-mode also applies)
#' @param maxShift x;      grid mode, horizontal distance to search over
#' @param step x;          grid mode, horizontal step size
#'
#  Optimiser mode operation (multi-mode also applies)
#' @param simplex         use simplex optimisation rather than doing the full bullseye plot
#' @param anneal          use simulated annealing optimisation
#' @param fixFsig         fix fSigma in simplex
#' @param geoError expError correlDist; rapid geolocation, using expected geolocation
#' error and correlation distance. Vertical shifts must be separatley defined
#' @param quickGeo        perform rapid geolocation using default error values
#' @param optTol x;        tolerance for optimisation
#' @param maxIter n;       maximum number of iterations
#' @param writeSimProg    write progress of simplex to output
#' @param writeWaves name; write out final waveforms as HDF5 when using simplex
#' @param nTriesAnneal n;  how many points do we try before stepping?
#' @param itersFixedT n;   how many iterations for each T?
#' @param kAnneal x;       Boltzmann constant for annealing
#' @param tInitial x; A initial annealing temperature
#' @param muAnneal x;      damping factor for temperature
#' @param tMinAnneal x;    minimum annealing temperature
#'
#  Initial estimates. Will search around this point
#' @param hOffset dx dy;   centre of horizontal offsets
#' @param offset z;        vertical datum offset
#'
#' Waveform characteristics
#' @param fSigma x;        footprint width, sigma in metres
#' @param pSigma x;        Gaussian pulse length, sigma in metres
#' @param readPulse file;  pulse shape, if not Gaussian
#'
#  Filters for input data
#' @param minSense x;      minimum LVIS/GEDI beam sensitivity to accept
#' @param maxZen zen;      maximum LVIS/GEDI zenith angle to use, degrees
#' @param maxScanAng ang;  maximum ALS scan angle, degrees
#' @param minDense x;      minimum ALS beam density to accept
#' @param decimate f;      decimate ALS point cloud by a factor, to save RAM
#' @param noFilt          don't filter outliers from correlation (default)
#' @param filtOut s;       filter outliers from correlation stats, along with
#' the number of standard deviations to use as a threshold
#' @param smooth sig;      smooth both waves before comparing
#' @param checkCover      only include footprints that are at least 2/3 covered with ALS data
#' @param median          use median correlation rather than mean
#'
#  Simulator settings. For simulator validation only
#' @param noNorm          don't correct sims for ALS densiy variations
#' @param norm            correct sims for ALS densiy variations
#' @param allSimMeth      use all simulation methods
#' @param pulseBefore     apply pulse shape before binning to prevent aliasing
#'
#  Octree to speed searching of ALS data. Not fully operational
#' @param noOctree       do not use an octree
#' @param octLevels n;    number of octree levels to use
#' @param nOctPix n;      number of octree pixels along a side for the top level
#'
#' 
#' @return Returns an S4 object of class [`hdf5r::H5File-class`]
#' containing the collocated waves
#'
#' @seealso
#' i) Hancock, S., Armston, J., Hofton, M., Sun, X., Tang, H., Duncanson, L.I., Kellner,
#' J.R. and Dubayah, R., 2019. The GEDI simulator: A large-footprint waveform lidar simulator
#' for calibration and validation of spaceborne missions. Earth and Space Science.
#' \doi{10.1029/2018EA000506}
#'
#' ii) gediSimulator: \url{https://bitbucket.org/StevenHancock/gedisimulator/src/master/}
#'
#' @examples
#' \dontshow{
#' rm(list = ls())
#' }
#'
#' @import hdf5r
#' @useDynLib rGEDI
#' @export
gediWFCollocate <- function(
    # Input-output
    output,
    als,
    gedi,
    aEPSG,
    readHDFgedi = TRUE,
    lgw = FALSE,
    readHDFlvis = FALSE,
    bounds = NULL,
    leaveEmpty = NULL,
    lEPSG = 4326,
    beamList = NULL,
    skipBeams = NULL,
    readBeams = NULL,
    # Multi-mode options
    solveCofG = NULL,
    maxVshift = -10.0,
    vStep = 0.2,
    # Grid mode operation (multi-mode also applies)
    maxShift = 10.0,
    step = 1.0,
    # Optimiser mode operation (multi-mode also applies)
    simplex = FALSE,
    anneal = FALSE,
    fixFsig = FALSE,
    geoError = FALSE,
    quickGeo = FALSE,
    optTol = 0.2,
    maxIter = 300,
    writeSimProg = NULL,
    writeWaves = NULL,
    nTriesAnneal = NULL,
    itersFixedT = NULL,
    kAnneal = NULL,
    tInitial = NULL,
    muAnneal = NULL,
    tMinAnneal = NULL,
    # Initial estimates. Will search around this point
    hOffset = NULL,
    offset = NULL,
    # Waveform characteristics
    fSigma = NULL,
    pSigma = NULL,
    readPulse = NULL,
    # Filters for input data
    minSense = NULL,
    maxZen = NULL,
    maxScanAng = NULL,
    minDense = NULL,
    decimate = NULL,
    noFilt = NULL,
    filtOut = NULL,
    smooth = NULL,
    checkCover = NULL,
    median = NULL,
    # Simulator settings. For simulator validation only
    noNorm = NULL,
    norm = NULL,
    allSimMeth = NULL,
    pulseBefore = NULL,
    # Octree to speed searching of ALS data. Not fully operational
    noOctree = NULL,
    octLevels = NULL,
    nOctPix = NULL) {
  # Check values
  stopifnotMessage(
    "als file(s) do not exist!" = all(file.exists(als)),
    "gedi file(s) do not exist!" = all(file.exists(gedi)),
    "output path is not valid!" = dir.exists(fs::path_dir(output)),
    "readPulse is invalid!" = checkFilepath(readPulse, newFile = FALSE, optional = TRUE),
    "fSigma is invalid!" = checkNumeric(fSigma),
    "checkCover is invalid!" = checkLogical(checkCover),
    "maxScanAng is invalid!" = checkNumeric(maxScanAng),
    "decimate is invalid!" = checkNumeric(decimate),
    "pulseBefore is invalid!" = checkLogical(pulseBefore),
    "noNorm is invalid!" = checkLogical(noNorm),
    "noOctree is invalid!" = checkLogical(noOctree),
    "octLevels is invalid!" = checkInteger(octLevels),
    "nOctPix is invalid!" = checkInteger(nOctPix)
  )

  gediInList <- inputOrInList(gedi)
  alsInList <- inputOrInList(als)
  if (fs::path_ext(output) != "h5") {
    output <- fs::path_ext_set(output, ".h5")
  }

  .Call(
    "C_collocateWaves",
    output,
    alsInList[[2]],
    alsInList[[1]],
    gediInList[[1]],
    gediInList[[2]],
    readHDFgedi,
    lgw,
    readHDFlvis,
    bounds,
    leaveEmpty,
    lEPSG,
    aEPSG,
    beamList,
    skipBeams,
    readBeams,

    # Multi-mode options
    solveCofG,
    maxVshift,
    vStep,

    # Grid mode operation (multi-mode also applies)
    maxShift,
    step,

    # Optimiser mode operation (multi-mode also applies)
    simplex,
    anneal,
    fixFsig,
    geoError,
    quickGeo,
    optTol,
    maxIter,
    writeSimProg,
    writeWaves,
    nTriesAnneal,
    itersFixedT,
    kAnneal,
    tInitial,
    muAnneal,
    tMinAnneal,

    # Initial estimates. Will search around this point
    hOffset,
    offset,

    # Waveform characteristics
    fSigma,
    pSigma,
    readPulse,

    # Filters for input data
    minSense,
    maxZen,
    maxScanAng,
    minDense,
    decimate,
    noFilt,
    filtOut,
    smooth,
    checkCover,
    median,

    # Simulator settings. For simulator validation only
    noNorm,
    norm,
    allSimMeth,
    pulseBefore,

    # Octree to speed searching of ALS data. Not fully operational
    noOctree,
    octLevels,
    nOctPix
  )

  unloadLibrary()

  cleanInList(alsInList)
  cleanInList(gediInList)
}

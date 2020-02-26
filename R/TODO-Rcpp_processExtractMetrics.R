#'GEDI full waveform data processing
#'
#'@description GEDI full waveform data processing and metrics extraction
#'
#'@param input waveform  input filename
# @param writeFit write fitted waveform
# @param writeGauss write Gaussian parameters
#' @param level2 level2 filename for LVIS ZG
#' @param bounds only analyse data within bounds.
#'
# Switches
# @param ground read true ground from file
# @param useInt use discrete intensity instead of count
# @param useFrac use fractional hits rather than counts
#' @param rhRes percentage energy resolution of RH metrics
#' @param laiRes LAI profile resolution in metres
#' @param laiH height to calculate LAI to
#' @param noRHgauss do not fit Gaussians
#' @param gTol ALS ground tolerance. Used to calculate slope.
#' @param fhdHistRes waveform intensity resolution to use when calculating FHD from histograms
#' @param forcePsigma do not read pulse sigma from file
#' @param bayesGround use Bayseian ground finding
# @param dontTrustGround don't trust ground in waveforms, if included
#' @param noRoundCoord do not round up coords when outputting

# Adding noise
#' @param dcBias mean noise level
#' @param nSig noise sigma
#' @param seed random number seed
#' @param hNoise hard threshold noise as a fraction of integral
#' @param linkNoise apply Gaussian noise based on link margin at a cover
#' @param linkFsig footprint width to use when calculating and applying signal noise
#' @param linkPsig pulse width to use when calculating and applying signal noise
#' @param trueSig true sigma of background noise
#' @param bitRate digitisation bit rate
#' @param maxDN maximum DN
#' @param renoise remove noise from truth before applying new noise level
#' @param newPsig new value for pulse width, when lengthening pulse
#' @param oldPsig old value for pulse width if not defined in waveform file, when lengthening pulse
#' @param addDrift apply detector background drift
#' @param missGround assume ground is missed to assess RH metrics
#' @param minGap delete signal beneath min detectable gap fraction

# Denoising
#' @param meanN mean noise level, if using a predefined mean level
#' @param thresh noise threshold, if using a predefined noise threshold
#' @param varNoise use a variable noise threshold
#' @param varScale variable noise threshold scale (multiple of stdev above mean to set threshold)
#' @param statsLen length to calculate noise stats over for varNoise
#' @param noiseTrack use noise tracking
#' @param sWidth smoothing width, after densoising
#' @param psWidth smoothing width, before denoising
#' @param msWidth smoothing width, after noise stats, before denoising
#' @param preMatchF matched filter before denoising
#' @param postMatchF matched filter after denoising
#' @param pFile read pulse file, for deconvolution and matched filters
#' @param gWidth Gaussian parameter selection smoothing width
#' @param minGsig minimum Gaussian sigma to fit
#' @param minWidth minimum feature width in bins
#' @param medNoise use median stats rather than mean
#' @param varDrift correct detector drift with variable factor
#' @param driftFac fix drift with constant drift factor
# @param rhoG ground reflectance
# @param rhoC canopy reflectance
#' @param pSigma pulse width to smooth by if using Gaussian pulse
#' @param gold deconvolve with Gold's method
#' @param deconTol deconvolution tolerance
#'
#' @examples
#'# GEDI level1B file path
#'level1b_filepath <- system.file("extdata", "gedi_level1b_clip.h5", package="rGEDI")
#'
#'FWmetrics =gediProcessFWaveform(level1b_filepath)
#'
# @useDynLib rGEDI
#' @import Rcpp methods
#' @export
gediProcessFWaveform = function(input, level2 = character(0), bounds = numeric(0), rhRes = 5, laiRes = 10.0, laiH = 30.0, noRHgauss = FALSE, gTol = 0.0, fhdHistRes = 0.001, forcePsigma = FALSE, bayesGround = FALSE, noRoundCoord = FALSE, dcBias = numeric(0), nSig = 0.0, seed = 1L, hNoise = 0.0, linkNoise = numeric(0), linkFsig = 5.5, linkPsig = 0.764331, trueSig = 5.0, bitRate = 12L, maxDN = 4096.0, renoise = FALSE, newPsig = 1.0, oldPsig = 0.764331, addDrift = 0.0, missGround = 0L, minGap = numeric(0), meanN = 0.0, thresh = 0.00000001, varNoise = TRUE, varScale = 3, statsLen = 10, noiseTrack = FALSE, sWidth = 0.3, psWidth = 0.0, msWidth = 0.0, preMatchF = FALSE, postMatchF = FALSE, pFile = character(0), gWidth = 1.2, minGsig = 0.764331, minWidth = 3, medNoise = FALSE, varDrift = FALSE, driftFac = numeric(0), pSigma = 1.0, gold = FALSE, deconTol = 0.0000001) {
  # Load dynamic library
  library.dynam(chname="rGEDI", package="rGEDI", lib.loc=.libPaths())

  if (h5::is.h5file(input)) {
    h5.file = h5::h5file(input, mode="r")
    datasets = h5::list.datasets(h5.file)
    expected = c(
      "/LON0",
      "/LAT0",
      "/LON1023",
      "/LAT1023",
      "/LFID",
      "/SHOTNUMBER",
      "/RXWAVE",
      "/INCIDENTANGLE",
      "/Z0",
      "/Z1023",
      "/SIGMEAN",
      "/TIME"
    )
    for (field in expected) {
      stopifnot(field %in% datasets)
    }
    h5::h5close(h5.file)
  } else {
    stop("input file does not seems to be valid HDF5 LVIS file")
  }
  df=processFloWave2(
    input = input,
    writeFit = FALSE,
    writeGauss = FALSE,
    readBinLVIS = FALSE,
    readHDFlvis = TRUE,
    readHDFgedi = FALSE,
    level2 = level2,
    bounds = bounds,
    ground = FALSE,
    useInt = FALSE,
    useFrac = FALSE,
    rhRes = rhRes,
    laiRes = laiRes,
    laiH = laiH,
    noRHgauss = noRHgauss,
    gTol = gTol,
    fhdHistRes = fhdHistRes,
    forcePsigma = forcePsigma,
    bayesGround = bayesGround,
    dontTrustGround = FALSE,
    noRoundCoord = noRoundCoord,
    dcBias = dcBias,
    nSig = nSig,
    seed = seed,
    hNoise = hNoise,
    linkNoise = linkNoise,
    linkFsig = linkFsig,
    linkPsig = linkPsig,
    trueSig = trueSig,
    bitRate = bitRate,
    maxDN = maxDN,
    renoise = renoise,
    newPsig = newPsig,
    oldPsig = oldPsig,
    addDrift = addDrift,
    missGround = missGround,
    minGap = minGap,
    meanN = meanN,
    thresh = thresh,
    varNoise = varNoise,
    varScale = varScale,
    statsLen = statsLen,
    noiseTrack = noiseTrack,
    sWidth = sWidth,
    psWidth = psWidth,
    msWidth = msWidth,
    preMatchF = preMatchF,
    postMatchF = postMatchF,
    pFile = pFile,
    gWidth = gWidth,
    minGsig = minGsig,
    minWidth = minWidth,
    medNoise = medNoise,
    varDrift = varDrift,
    driftFac = driftFac,
    pSigma = pSigma,
    gold = gold,
    deconTol = deconTol
  )

  # Unload dynamic library
  library.dynam.unload("rGEDI", libpath=system.file(package="rGEDI"))
  return(df)
}

#'GEDI full waveform data processing
#'
#' @description GEDI full waveform data processing and metrics extraction
#'
#! \bold{Input output}
#' @param input name. waveform  input filename
#' @param outRoot name. output filename root
# @param inList list. input file list for multiple files
#' @param writeFit write fitted waveform
#' @param writeGauss write Gaussian parameters

# Always GEDI
# @param readBinLVIS input is an LVIS binary file
# @param readHDFlvis read LVIS HDF5 input
# @param readHDFgedi read GEDI simulator HDF5 input
# @param level2 name. level2 filename for LVIS ZG
#' @param bounds minX minY maxX maxY. only analyse data within bounds

# Not in use yet
# @param beamList character. 0/1 for whether or not to use beams 18, default "11111111"
# @param skipBeams character. list of beam numbers to skip. No spaces between (eg "123")
# @param readBeams character. list of beam numbers to read. No spaces between (eg "123")
#'
#! \bold{Switches}

# Always needed
# @param ground read true ground from file
#' @param useInt use discrete intensity instead of count
#' @param useFrac use fractional hits rather than counts
#' @param rhRes r. percentage energy resolution of RH metrics
#' @param laiRes res. lai profile resolution in metres
#' @param laiH h. height to calculate LAI to
#' @param noRHgauss do not fit Gaussians
#' @param gTol tol. ALS ground tolerance. Used to calculate slope.
#' @param fhdHistRes res. waveform intensity resolution to use when calculating FHD from histograms
#' @param forcePsigma do not read pulse sigma from file
#' @param bayesGround use Bayseian ground finding
#' @param dontTrustGround don't trust ground in waveforms, if included
#' @param noRoundCoord do not round up coords when outputting
#' @param noCanopy do not calculate FHD histograms and LAI profiles
#'
#! \bold{Adding noise}
#' @param dcBias n. mean noise level
#' @param nSig sig. noise sigma
#' @param seed n. random number seed
#' @param hNoise n. hard threshold noise as a fraction of integral
#' @param linkNoise linkM cov. apply Gaussian noise based on link margin at a cover
#' @param linkFsig sig. footprint width to use when calculating and applying signal noise
#' @param linkPsig sig. pulse width to use when calculating and applying signal noise
#' @param trueSig sig. true sigma of background noise
#' @param bitRate n. digitisation bit rate
#' @param maxDN max. maximum DN
#' @param renoise remove noise from truth before applying new noise level
#' @param newPsig sig. new value for pulse width, when lengthening pulse
#' @param oldPsig sig. old value for pulse width if not defined in waveform file, when lengthening pulse
#' @param addDrift xi. apply detector background drift
#' @param missGround assume ground is missed to assess RH metrics
#' @param minGap gap. delete signal beneath min detectable gap fraction
#'
#! \bold{Photon counting}
#' @param photonCount output point cloud from photon counting
#' @param pcl convert to photon counting pulsecompressed
#' @param nPhotons n. mean number of photons
#' @param photonWind x. window length for photon counting search, metres
#' @param noiseMult x. noise multiplier for photoncounting
#' @param rhoVrhoG x. ratio of canopy to ground reflectance at this wavelength. Not different from rhoV and rhoG
#' @param nPhotC n. mean number of canopy photons (replaces nPhotons and rhoVrhoG)
#' @param nPhotG n. mean number of ground photons (replaces nPhotons and rhoVrhoG)
#' @param photHDF write photoncounting
#'
#! \bold{Denoising}
#' @param meanN n. mean noise level, if using a predefined mean level
#' @param thresh n. noise threshold, if using a predefined noise threshold
#' @param varNoise use a variable noise threshold
#' @param varScale x. variable noise threshold scale (multiple of stdev above mean to set threshold)
#' @param statsLen len. length to calculate noise stats over for varNoise
#' @param noiseTrack use noise tracking
#' @param sWidth sig. smoothing width, after densoising
#' @param psWidth sigma. smoothing width, before denoising
#' @param msWidth sig. smoothing width, after noise stats, before denoising
#' @param preMatchF matched filter before denoising
#' @param postMatchF matched filter after denoising
#' @param pFile file. read pulse file, for deconvoltuion and matched filters
#' @param gWidth sig. Gaussian parameter selection smoothing width
#' @param minGsig sig. minimum Gaussian sigma to fit
#' @param minWidth n. minimum feature width in bins
#' @param medNoise use median stats rather than mean
#' @param varDrift correct detector drift with variable factor
#' @param driftFac xi. fix drift with constant drift factor
#' @param rhoG rho. ground reflectance
#' @param rhoC rho. canopy reflectance
#' @param pSigma sig. pulse width to smooth by if using Gaussian pulse
#' @param gold deconvolve with Gold's method
#' @param deconTol deconvolution tolerance
#'
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
#'# Simulate GEDI full-waveform
#'wf<-gediWFSimulator(input=LASfile,output="gediSimulation.h5")
#'
#'# Adding noise to GEDI full-waveform
#'wfn<-gediWFNoise(input=wf,output="gediSimulation_noise.h5")
#'
#'# Extracting GEDI feull-waveform derived metrics
#'wfmetrics<-gediWFMetrics(input=wfn,outRoot=getwd())
#'
#' @useDynLib rGEDI
#' @import methods
#' @export
gediWFMetric = function(
  input,
  outRoot,
  writeFit = FALSE,
  writeGauss = FALSE,
  bounds = NULL,
  # beamList = NULL,
  # skipBeams = NULL,
  # readBeams = NULL,
  ground = FALSE,
  useInt = FALSE,
  useFrac = FALSE,
  rhRes = 5.0,
  laiRes = 10.0,
  laiH = 30.0,
  noRHgauss = FALSE,
  gTol = 0.0,
  fhdHistRes = 0.001,
  forcePsigma = FALSE,
  bayesGround = FALSE,
  dontTrustGround = FALSE,
  noRoundCoord = FALSE,
  noCanopy = FALSE,
  dcBias = 0.0,
  nSig = 0.0,
  seed = NULL,
  hNoise = 0.0,
  linkNoise = NULL,
  linkFsig = NULL,
  linkPsig = NULL,
  trueSig = NULL,
  bitRate = NULL,
  maxDN = NULL,
  renoise = FALSE,
  newPsig = -1.0,
  oldPsig = 0.764331,
  addDrift = NULL,
  missGround = FALSE,
  minGap = NULL,
  photonCount = FALSE,
  pcl = FALSE,
  nPhotons = 2.1,
  photonWind = 200.0,
  noiseMult = 0.1,
  rhoVrhoG = 1.0,
  nPhotC = 2.1,
  nPhotG = -1.0,
  photHDF = FALSE,
  meanN = 0.0,
  thresh = 0.00000000000001,
  varNoise = FALSE,
  varScale = NULL,
  statsLen = NULL,
  noiseTrack = FALSE,
  sWidth = NULL,
  psWidth = 0.0,
  msWidth = NULL,
  preMatchF = FALSE,
  postMatchF = FALSE,
  pFile = NULL,
  gWidth = 1.2,
  minGsig = 0.764331,
  minWidth = 0.0,
  medNoise = FALSE,
  varDrift = NULL,
  driftFac = NULL,
  rhoG = 0.4,
  rhoC = 0.57,
  pSigma = NULL,
  gold = FALSE,
  deconTol = NULL) {

  readBinLVIS = FALSE
  readHDFlvis = FALSE
  readHDFgedi = TRUE
  level2 = NULL
  beamList = NULL
  skipBeams = NULL
  readBeams = NULL
  ground = FALSE

  stopifnotMessage(
    all(file.exists(input)),
    all(fs::path_ext(input) == "h5"),
    checkParentDir(outRoot, optional=FALSE),
    checkLogical(writeFit),
    checkLogical(writeGauss),
    checkNumericLength(bounds, 4),
    checkLogical(useInt),
    checkLogical(useFrac),
    checkNumeric(laiRes),
    checkNumeric(laiH),
    checkLogical(noRHgauss),
    checkNumeric(gTol),
    checkNumeric(fhdHistRes),
    checkLogical(forcePsigma),
    checkLogical(bayesGround),
    checkLogical(dontTrustGround),
    checkLogical(noRoundCoord),
    checkLogical(noCanopy),
    checkNumeric(dcBias),
    checkNumeric(nSig),
    checkInteger(seed),
    checkNumeric(hNoise),
    checkNumericLength(linkNoise, 2),
    checkNumeric(linkFsig),
    checkNumeric(linkPsig),
    checkNumeric(trueSig),
    checkInteger(bitRate),
    checkNumeric(maxDN),
    checkLogical(renoise),
    checkNumeric(newPsig),
    checkNumeric(oldPsig),
    checkNumeric(addDrift),
    checkLogical(missGround),
    checkLogical(minGap),
    checkLogical(photonCount),
    checkLogical(pcl),
    checkNumeric(nPhotons),
    checkNumeric(photonWind),
    checkNumeric(noiseMult),
    checkNumeric(rhoVrhoG),
    checkNumeric(nPhotC),
    checkNumeric(nPhotG),
    checkLogical(photHDF),
    checkNumeric(meanN),
    checkNumeric(thresh),
    checkLogical(varNoise),
    checkNumeric(varScale),
    checkNumeric(statsLen),
    checkLogical(noiseTrack),
    checkNumeric(sWidth),
    checkNumeric(psWidth),
    checkNumeric(msWidth),
    checkLogical(preMatchF),
    checkLogical(postMatchF),
    checkFilepath(pFile, newFile=FALSE, optional=TRUE),
    checkNumeric(gWidth),
    checkNumeric(minGsig),
    checkNumeric(minWidth),
    checkLogical(medNoise),
    checkLogical(varDrift),
    checkNumeric(driftFac),
    checkNumeric(rhoG),
    checkNumeric(rhoC),
    checkNumeric(pSigma),
    checkLogical(gold),
    checkNumeric(deconTol)
  )

  inputInList = inputOrInList(input)

  res = .Call("C_gediMetrics",
        # Input output
        inputInList[[1]],
        outRoot,
        inputInList[[2]],
        writeFit,
        writeGauss,
        readBinLVIS,
        readHDFlvis,
        readHDFgedi,
        level2,
        bounds,
        beamList,
        skipBeams,
        readBeams,

        # Switches
        ground,
        useInt,
        useFrac,
        rhRes,
        laiRes,
        laiH,
        noRHgauss,
        gTol,
        fhdHistRes,
        forcePsigma,
        bayesGround,
        dontTrustGround,
        noRoundCoord,
        noCanopy,

        # Adding noise
        dcBias,
        nSig,
        seed,
        hNoise,
        linkNoise,
        linkFsig,
        linkPsig,
        trueSig,
        bitRate,
        maxDN,
        renoise,
        newPsig,
        oldPsig,
        addDrift,
        missGround,
        minGap,

        # Photon Counting
        photonCount,
        pcl,
        nPhotons,
        photonWind,
        noiseMult,
        rhoVrhoG,
        nPhotC,
        nPhotG,
        photHDF,

        # Denoising
        meanN,
        thresh,
        varNoise,
        varScale,
        statsLen,
        noiseTrack,
        sWidth,
        psWidth,
        msWidth,
        preMatchF,
        postMatchF,
        pFile,
        list(gWidth,
             minGsig,
             minWidth,
             medNoise,
             varDrift,
             driftFac,
             rhoG,
             rhoC,
             pSigma,
             gold,
             deconTol))
  unloadLibrary()
  cleanInList(inputInList)

  if (res==0) {
    output = paste0(outRoot, ".metric.txt")
    header = read.csv(output, sep=",", nrow=1, header = FALSE, as.is=TRUE)
    header[ncol(header)] = NULL
    header=gsub(" *\\d+ *([^ ]*)", "\\1", header)
    metricData = read.csv(output, sep=" ", skip=1, na.strings = "?", header = FALSE)
    if (ncol(metricData) == length(header)) {
      names(metricData) = header
    }
  }

  return (metricData)
}
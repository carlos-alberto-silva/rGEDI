# In R ----------------------------------------
#' @useDynLib rGEDI
#' @export
gediSimulator <- function(
  input,
  output,
  inList = NULL,
  ground = FALSE,
  hdf = TRUE,
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
  LVIS = FALSE,
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

  decon = FALSE,
  indDecon = FALSE,
  readWave = FALSE,

  listFiles = FALSE,
  keepOld = FALSE,
  useShadow = FALSE,
  polyGround = FALSE,
  nnGround = FALSE,
  seed = NULL) {
  .Call("C_gedisimulate",
      input,
      output,
      inList,
      ground,
      hdf,
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
      seed)
  return(hdf5r::H5File$new(output, "r+"))
}
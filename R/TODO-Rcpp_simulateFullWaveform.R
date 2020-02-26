#' GEDI full waveform data simulation
#'
#' @description Simulate GEDI full waveform data from Airborne Laser Scanning (ALS) 3-D point cloud
#'
#  Input output filenames and format
#' @param input name. lasfile input filename
#' @param output name. output filename
#' @param coords lon lat matrix. footprint coordinate in same system as lasfile
#' @param inList list. input file list (ASCII file) for multiple files
#' @param ground record separate ground and canopy waveforms
#' @param hdf write output as HDF5, default TRUE. Best with gridded or list of coords
#' @param waveID id.
#
#
#  Lidar characteristics. Defaults are expected GEDI values.
#' @param pSigma sig. set Gaussian pulse width as 1 sigma
#' @param pFWHM fhwm. set Gaussian pulse width as FWHM in ns
#' @param readPulse file. read pulse shape and width from a file instead of making Gaussian
#' @param fSigma sig. set footprint width
#' @param wavefront file. read wavefront shape from file instead of setting Gaussian. Note that footprint width is still set by fSigma
#' @param res res. range resolution of waveform digitisation to output, in units of ALS data
#' @param LVIS use LVIS pulse length, sigma=6.25m
#' @param topHat use a top hat wavefront
#' @param sideLobe use side lobes
#' @param lobeAng ang. lobe axis azimuth
#
#  Input data quality filters
#' @param checkCover check that the footprint is covered by ALS data. Do not output if not
#' @param maxScanAng ang. maximum scan angle, degrees
#' @param decimate x. probability of accepting an ALS beam

#  Computational speed options
#' @param pBuff s. point reading buffer size in Gbytes
#' @param maxBins Optional: for HDF5, limit number of bins to save trimming.
#' @param countOnly only use count method
#' @param pulseAfter apply the pulse smoothing after binning for computational speed, at the risk of aliasing (default)
#' @param pulseBefore apply the pulse smoothing before binning to avoid the risk of aliasing, at the expense of computational speed
#' @param noNorm don't normalise for ALS density
#
#  Octree
#' @param noOctree do not use an octree
#' @param octLevels n. number of octree levels to use
#' @param nOctPix n. number of octree pixels along a side for the top level
#
#  Using full-waveform input data (not tested)
#' @param decon deconvolve
#' @param indDecon deconvolve individual beams
#' @param readWave read full-waveform where available
#
#  Miscellaneous
#' @param listFiles list files. Do not read them
#' @param keepOld do not overwrite old files, if they exist
#' @param useShadow account for shadowing in discrete return data through voxelisation
#' @param polyGround find mean ground elevation and slope through fitting a polynomial
#' @param nnGround find mean ground elevation and slope through nearest neighbour
#' @param seed n. seed number for random numbers generator
#'
# @useDynLib rGEDI
#' @import Rcpp methods h5
#' @export
gediFWSimulation = function(input, output, coords, inList = FALSE, ground = FALSE, hdf = TRUE, waveID = character(0), pSigma = -1.0, pFWHM = 15.0, readPulse = character(0), fSigma = 5.5, wavefront = character(0), res = 0.15, LVIS = FALSE, topHat = FALSE, sideLobe = FALSE, lobeAng = 0.0, checkCover = FALSE, maxScanAng = 1000000.0, decimate = 1.0, pBuff = 0.2, maxBins = 1024L, countOnly = FALSE, pulseAfter = FALSE, pulseBefore = FALSE, noNorm = FALSE, noOctree = FALSE, octLevels = 0L, nOctPix = 40L, decon = FALSE, indDecon = FALSE, readWave = FALSE, listFiles = FALSE, keepOld = FALSE, useShadow = FALSE, polyGround = FALSE, nnGround = FALSE, seed = numeric(0)) {
  # Load dynamic library
  library.dynam(chname="rGEDI", package="rGEDI", lib.loc=.libPaths())
  coordsList=tempfile()
  write.table(coords, file = coordsList, row.names = FALSE, col.names = FALSE)
  gediRat(
    input,
    output,
    inList,
    ground,
    hdf=hdf,
    waveID=waveID,
    coords=numeric(0),
    listCoord=coordsList,
    gridBound=numeric(0),
    gridStep = 30,
    pSigma=pSigma,
    pFWHM=pFWHM,
    readPulse=readPulse,
    fSigma=fSigma,
    wavefront=wavefront,
    res=res,
    LVIS=LVIS,
    topHat=topHat,
    sideLobe=sideLobe,
    lobeAng=lobeAng,
    checkCover=checkCover,
    maxScanAng=maxScanAng,
    decimate=decimate,
    pBuff=pBuff,
    maxBins=maxBins,
    countOnly=countOnly,
    pulseAfter=pulseAfter,
    pulseBefore=pulseBefore,
    noNorm=noNorm,
    noOctree=noOctree,
    octLevels=octLevels,
    nOctPix=nOctPix,
    decon=decon,
    indDecon=indDecon,
    readWave=readWave,
    listFiles=listFiles,
    keepOld=keepOld,
    useShadow=useShadow,
    polyGround=polyGround,
    nnGround=nnGround,
    seed = numeric(0)
  )

  # Unload dynamic library
  library.dynam.unload("rGEDI", libpath=system.file(package="rGEDI"))

  message(paste0("Successfully written to ", output))
}

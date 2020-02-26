#' Add noise to simulated GEDI waveform
#'
#' @description
#' This function is intended to be used after the simulation process (\code{\link{gediWFSimulator}}).
#' Given the known laser power, optical efficiences, mean atmospheric
#' transmission at 1,064 nm, expected canopy and ground reflectance,
#' range of background illumination intensities, and the detector response
#' the expected performance is calculated according to \href{https://doi.org/10.1109/26.8924}{Davidson and Sun (1988)}.
#'
#' @param input \code{\link{character}}. Waveform input filename
#' @param output \code{\link{character}}. Output filename
#' @param seed \code{\link{integer}}. Random number seed
#' @param linkNoise \code{\link{numeric}(2)}. Apply Gaussian noise based on link margin at a cover
#' @param dcBias \code{\link{numeric}}. Mean noise level
#' @param linkFsig \code{\link{numeric}}. Footprint width to use when calculating and applying signal noise
#' @param linkPsig \code{\link{numeric}}. Pulse width to use when calculating and applying signal noise
#' @param trueSig \code{\link{numeric}}. True sigma of background noise
#' @param bitRate \code{\link{integer}}. DN bit rate
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
#'#'# specify the path to ALS data
#'LASfile <- system.file("extdata", "LASexample1.las", package="rGEDI")
#'
#'# Simulate GEDI full-waveform
#'wf<-gediWFSimulator(input=LASfile,output="gediSimulation.h5")
#'
#'# Adding noise to GEDI full-waveform
#'wfn<-gediWFNoise(input=wf,output="gediSimulation_noise.h5")
#'
#' @export
gediWFNoise <- function(
  input,
  output,
  seed = NULL,
  linkNoise = c(3.0103, 0.95),
  dcBias = NULL,
  linkFsig = NULL,
  linkPsig = 6.6383,
  trueSig = NULL,
  bitRate = 12) {



  # Check values
  stopifnotMessage(
    all(file.exists(input)),
    all(fs::path_ext(input) == "h5"),
    checkFilepath(output, newFile=TRUE, optional=FALSE),
    checkInteger(seed),
    checkNumericLength(linkNoise, 2),
    checkNumeric(dcBias),
    checkNumeric(linkFsig),
    checkNumeric(linkPsig),
    checkNumeric(trueSig),
    checkInteger(bitRate)
  )
  res = .Call("C_addNoiseHDF",
        input,
        output,
        seed,
        linkNoise,
        dcBias,
        linkFsig,
        linkPsig,
        trueSig,
        bitRate)
  unloadLibrary()
  return(hdf5r::H5File$new(output, "r+"))
}

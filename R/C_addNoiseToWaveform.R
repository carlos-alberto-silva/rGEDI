#' Add noise to simulated GEDI full-waveform
#'
#' @description This function adds noise to the simulated GEDI full-waveform (output of gediWFSimulator)
#'
#' @usage gediWFNoise(input,output,seed,linkNoise,dcBias,linkFsig,linkPsig,trueSig,bitRate)
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
  linkNoise = NULL,
  dcBias = NULL,
  linkFsig = NULL,
  linkPsig = NULL,
  trueSig = NULL,
  bitRate = NULL) {
  loadLibrary()
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

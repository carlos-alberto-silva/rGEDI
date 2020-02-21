#' GEDI add noise to simulated waveform
#'
#' @description Add noise to HDF output from gediWFSimulator
#'
#' @param input \code{\link{character}}. Waveform  input filename
#' @param output \code{\link{character}}. Output filename
#' @param seed \code{\link{integer}}. Random number seed
#' @param linkNoise \code{\link{numeric}(2)}. Apply Gaussian noise based on link margin at a cover
#' @param dcBias \code{\link{numeric}}. Mean noise level
#' @param linkFsig \code{\link{numeric}}. Footprint width to use when calculating and applying signal noise
#' @param linkPsig \code{\link{numeric}}. Pulse width to use when calculating and applying signal noise
#' @param trueSig \code{\link{numeric}}. True sigma of background noise
#' @param bitRate \code{\link{integer}}. DN bit rate
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

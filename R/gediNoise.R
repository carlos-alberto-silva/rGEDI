#' @export
gediNoise <- function(
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

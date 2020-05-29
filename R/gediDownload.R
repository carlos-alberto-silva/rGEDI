#'Download GEDI data
#'
#'@description Download GEDI data from LP DAAC Data Pool. Users will need to enter their
#'Earth Explore login Information for downloading the data.
#'
#'@param filepath Vector object; path to the GEDI data
#'@param outdir Vector object, output directory for downloading GEDI data, default tempdir()
#'@param overwrite logical; overwrite file if they already exists in destination, default FALSE
#'@param buffer_size integer; the size of download chunk in KB to hold in memory before writing to file, default 512.
#' 
#'@return No return value on success, on failure it will \code{stop()}
#'@references Credits to Cole Krehbiel. Code adapted from <https://git.earthdata.nasa.gov/projects/LPDUR/repos/daac_data_download_r/browse/DAACDataDownload.R>
#'@examples
#'\donttest{
#'# Set path to GEDI data
#'# herein we will only download xml metedata
#'filepath=c(paste0(
#'                  "https://e4ftl01.cr.usgs.gov/GEDI/GEDI02_B.001",
#'                  "/2019.04.18/GEDI02_B_2019108032534_O01961_T03911_02_001_01.h5.xml"
#'                  ),
#'           paste0("https://e4ftl01.cr.usgs.gov/GEDI/GEDI02_B.001",
#'                  "/2019.04.18/GEDI02_B_2019108045815_O01962_T01066_02_001_01.h5.xml"
#'                 )
#'          )
#'
#'# Set dir to download files to
#'outdir=tempdir()
#'
#'# Create .netrc file
#'netrc = file.path(outdir, ".netrc")
#'netrc_conn <- file(netrc)
#'
#'writeLines(c("machine urs.earthdata.nasa.gov",
#'             sprintf("login %s", Sys.getenv("NASA_USER")),
#'             sprintf("password %s", Sys.getenv("NASA_PASSWORD"))
#'), netrc_conn)
#'
#'close(netrc_conn)
#'
#'#' Downloading GEDI data
#'gediDownload(filepath,outdir)
#'}
#'@import curl
#'@export
gediDownload<-function(filepath, outdir = NULL, overwrite = FALSE, buffer_size = 512){
  if (is.null(outdir)) {
    outdir == tempdir()
  }
  stopifnotMessage(
    "outdir is not a valid path" = checkParentDir(outdir),
    "overwrite is not logical" = checkLogical(overwrite),
    "buffer_size is not an integer" = checkInteger(buffer_size)
  )
  buffer_size = as.integer(buffer_size)
  netrc = getNetRC(outdir)

  files<-filepath
  n_files = length(files)

  # Download all files in filepath vector
  for (i in 1:n_files) {
    url = files[i]
    message("------------------------------")
    message(sprintf("Downloading file %d/%d: %s", i, n_files, basename(url)))
    message("------------------------------")

    if (gediDownloadFile(
      url,
      outdir,
      overwrite,
      buffer_size,
      netrc
    ) == 0) {
      message("Finished successfully!")
    } else {
      stop(sprintf("File %s has not been downloaded properly!", basename(url)))
    }
  }
}

gediDownloadFile = function(url, outdir, overwrite, buffer_size, netrc) {
  filename <- file.path(outdir, basename(url)) # Keep original filename
  if((! overwrite) && file.exists(filename)) {
    message("Skipping this file, already downloaded!")
    return(0)
  } # SKip if already downloaded

  # Temporary to file to resume to
  resume=paste0(filename, ".curltmp")
  if(file.exists(resume)){
    resume_from=file.info(resume)$size # Get current size to resume from
  } else {
    resume_from=0
  }

  # Connection config
  h = curl::new_handle()
  curl::handle_setopt(h, netrc=1, netrc_file=netrc, resume_from=resume_from)

  tryCatch({
    fileHandle=file(resume, open="ab", raw = T)
    conn=curl::curl(url, handle=h, open="rb")
    headers=rawToChar(curl::handle_data(h)$headers)
    total_size=as.numeric(gsub("[^\u00e7]*Content-Length: ([0-9]+)[^\u00e7]*","\\1",x=headers, perl = T))
    while(TRUE) {
      message(sprintf("\rDownloading... %.2f/%.2fMB (%.2f%%)    ",
                      resume_from/1024.0/1024.0,
                      total_size/1024.0/1024.0,
                      100.0*resume_from/total_size),
              appendLF=FALSE)
      data = readBin(conn, what = raw(), n = 1024*buffer_size)
      size = length(data)
      if (size==0) {
        break
      }
      writeBin(data, fileHandle, useBytes = T)
      resume_from = resume_from + size
    }
    message(sprintf("\rDownloading... %.2f/%.2fMB (100%%)    ",
            total_size/1024.0/1024.0,
            total_size/1024.0/1024.0))
    close(fileHandle)
    close(conn)
    file.rename(resume, filename)
    return(0)
  }, interrupt=function(e){
    warning("\nDownload interrupted!!!")
    try(close(conn), silent = TRUE)
    try(close(fileHandle), silent = TRUE)
  }, finally = {
    try(close(conn), silent = TRUE)
    try(close(fileHandle), silent = TRUE)
  })
  return(-1)
}

getNetRC = function(dl_dir) {
  netrc <- file.path(dl_dir,'.netrc')  # Path to netrc file
  # ------------------------------------CREATE .NETRC FILE------------------------------------------ #
  if (file.exists(netrc) == FALSE || grepl("urs.earthdata.nasa.gov", readLines(netrc)) == FALSE) {
    netrc_conn <- file(netrc)

    # User will be prompted for NASA Earthdata Login Username and Password below
    writeLines(c("machine urs.earthdata.nasa.gov",
                 sprintf("login %s", getPass::getPass(msg = "Enter NASA Earthdata Login Username \n (or create an account at urs.earthdata.nasa.gov) :")),
                 sprintf("password %s", getPass::getPass(msg = "Enter NASA Earthdata Login Password:"))), netrc_conn)
    close(netrc_conn)
    message("A .netrc file with your Earthdata Login credentials was stored in the output directory ")
  }
  return (netrc)
}

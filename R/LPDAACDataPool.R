#'Download GEDI data
#'
#'@description Download GEDI data from LP DAAC Data Pool. Users will need to enter their
#'Earth Explore login Information for downloading the data.
#'
#'@usage LPDAACDataPool(filepath,outdir)
#'
#'@param filepath Vector object; path to the GEDI data
#'@param outdir Vector object, output directory for downloading GEDI data
#'@references Credits to Cole Krehbiel. Code adpted from <https://git.earthdata.nasa.gov/projects/LPDUR/repos/daac_data_download_r/browse/DAACDataDownload.R>
#'@examples
#'\donttest{
#'#' Set path to GEDI data
#'filepath=c(paste0(
#'                  "https://e4ftl01.cr.usgs.gov/GEDI/GEDI02_B.001",
#'                  "/2019.04.18/GEDI02_B_2019108032534_O01961_T03911_02_001_01.h5"
#'                  ),
#'           paste0("https://e4ftl01.cr.usgs.gov/GEDI/GEDI02_B.001",
#'                  "/2019.04.18/GEDI02_B_2019108045815_O01962_T01066_02_001_01.h5"
#'                 )
#'          )
#'
#'# Set dir to download files to
#'outdir=tempdir()
#'
#'#' Downloading GEDI data
#'LPDAACDataPool(filepath,outdir)
#'}
#'@export
LPDAACDataPool<-function(filepath,outdir){
    # ---------------------------------SET UP ENVIRONMENT--------------------------------------------- #
    # IMPORTANT: Update the line below if you want to download to a different directory (ex: "c:/data/")
    if (is.null(outdir)){
      dl_dir <- tempdir()                                           # Set dir to download files to
    } else {
      dl_dir <- outdir                                           # Set dir to download files to
    }
    netrc <- file.path(dl_dir,'.netrc', fsep = .Platform$file.sep)  # Path to netrc file
    # ------------------------------------CREATE .NETRC FILE------------------------------------------ #
    if (file.exists(netrc) == FALSE || grepl("urs.earthdata.nasa.gov", readLines(netrc)) == FALSE) {
      netrc_conn <- file(netrc)

      # User will be prompted for NASA Earthdata Login Username and Password below
      writeLines(c("machine urs.earthdata.nasa.gov",
                   sprintf("login %s", getPass::getPass(msg = "Enter NASA Earthdata Login Username \n (or create an account at urs.earthdata.nasa.gov) :")),
                   sprintf("password %s", getPass::getPass(msg = "Enter NASA Earthdata Login Password:"))), netrc_conn)
      close(netrc_conn)
      cat("A .netrc file with your Earthdata Login credentials was stored in the output directory \n")
    }

    # ---------------------------CONNECT TO DATA POOL AND DOWNLOAD FILES------------------------------ #
    files<-filepath
    # Loop through all files
    for (i in 1:length(files)) {
      filename <-  tail(strsplit(files[i], '/')[[1]], n = 1) # Keep original filename

      cat(paste0("Downloading file ",filename," \n"))
      # Write file to disk (authenticating with netrc) using the current directory/filename
      response <- httr::GET(files[i], httr::write_disk(paste0(outdir,"\\",filename), overwrite = TRUE), httr::progress(),
                            httr::config(netrc = TRUE, netrc_file = netrc), httr::set_cookies("LC" = "cookies"))

      # Check to see if file downloaded correctly
      if (response$status_code == 200) {
        cat(sprintf("Done!","%s downloaded at %s", filename, dl_dir))
      } else {
        cat(sprintf("%s not downloaded. Verify that your username and password are correct in %s", filename, netrc))
      }
     }
}

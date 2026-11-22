#'Get GEDI Pulse Full Waveform (GEDI Level1B)
#'
#'@description This function extracts the full waveform of a given pulse from GEDI Level1B data.
#'
#'@usage getLevel1BWF(level1b, shot_number)
#'
#'@param level1b A GEDI Level1B object (output of [readLevel1B()] function). A S4 object of class "gedi.level1b".
#'@param shot_number Shot number. A scalar representing the shot number of a giving pulse.
#'
#'@return Returns an S4 object of class "gedi.fullwaveform".
#'
#'@details Shot numbers can be extracted using [readLevel1B] function.
#'
#'@seealso \url{https://www.earthdata.nasa.gov/data/catalog/lpcloud-gedi01-b-002}
#'
#'@examples
#'# Specifying the path to GEDI level1B data (zip file)
#'outdir = tempdir()
#'level1B_fp_zip <- system.file("extdata",
#'                   "GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level1B data
#'level1Bpath <- unzip(level1B_fp_zip,exdir = outdir)
#'
#'# Reading GEDI level1B data (h5 file)
#'level1b<-readLevel1B(level1Bpath=level1Bpath)
#'
#'# Extracting GEDI full waveform for a giving shotnumber
#'wf <- getLevel1BWF(level1b, shot_number="19640521100108408")
#'
#'# Plotting GEDI Full waveform
#'oldpar<-par()
#'par(mfrow = c(1,2), cex.axis = 1.5)
#'plot(wf, relative=FALSE, polygon=TRUE, type="l", lwd=2, col="forestgreen",
#'xlab="Waveform Amplitude", ylab="Elevation (m)")
#'
#'plot(wf, relative=TRUE, polygon=TRUE, type="l", lwd=2, col="forestgreen",
#'xlab="Waveform Amplitude (%)", ylab="Elevation (m)")
#'
#'par(oldpar)
#'close(level1b)
#'@export
#'
getLevel1BWF<-function(level1b,shot_number){

  if (inherits(level1b, "gedi_granule_collection")) {
    for (source in level1b$sources) {
      object <- readLevel1B(source)
      result <- tryCatch(getLevel1BWF(object, shot_number), error = identity)
      close(object)
      if (!inherits(result, "error")) return(result)
    }
    stop("Shot number ", shot_number, " was not found in the granule collection.",
         call. = FALSE)
  }

  level1b<-level1b@h5
  groups_id<-grep("BEAM\\d{4}$",gsub("/","",
                                     .gedi_list_groups(level1b, recursive = FALSE)), value = T)

  shot_text <- as.character(shot_number)
  beam_id <- suppressWarnings(as.integer(substr(shot_text, nchar(shot_text)-12L, nchar(shot_text)-11L)))
  if (length(beam_id) == 1L && !is.na(beam_id) && beam_id >= 0L && beam_id <= 15L) {
    preferred <- paste0("BEAM", paste(rev(as.integer(intToBits(beam_id))[1:4]), collapse = ""))
    groups_id <- c(intersect(preferred, groups_id), setdiff(groups_id, preferred))
  }
  i=NULL
  for ( k in groups_id){
    gid<-max(level1b[[paste0(k,"/shot_number")]][]==shot_number)
    if (gid==1) {i=k; break}
  }

  if(is.null(i)) {
    stop(paste0("Shot number ", shot_number, " was not found within the dataset!. Please try another shot number"))
  } else{

  shot_number_i<-level1b[[paste0(i,"/shot_number")]][]
  shot_number_id<-which(shot_number_i[]==shot_number)
  elevation_bin0<-level1b[[paste0(i,"/geolocation/elevation_bin0")]][]
  elevation_lastbin<-level1b[[paste0(i,"/geolocation/elevation_lastbin")]][]
  rx_sample_count<-level1b[[paste0(i,"/rx_sample_count")]][]
  # HDF5 addresses fit in an R double; integer64 cannot be used with colon.
  rx_sample_start_index<-as.numeric(level1b[[paste0(i,"/rx_sample_start_index")]][])
  rx_sample_start_index_n<-rx_sample_start_index-min(rx_sample_start_index)+1
  rxwaveform_i<-level1b[[paste0(i,"/rxwaveform")]][rx_sample_start_index_n[shot_number_id]:(rx_sample_start_index_n[shot_number_id]+rx_sample_count[shot_number_id]-1)]
  if (length(rxwaveform_i) != as.numeric(rx_sample_count[shot_number_id])) {
    stop("The requested waveform sample range is incomplete.", call. = FALSE)
  }
  rxwaveform_inorm<-(rxwaveform_i-min(rxwaveform_i))/(max(rxwaveform_i)-min(rxwaveform_i))*100
  elevation_bin0_i<-elevation_bin0[shot_number_id]
  elevation_lastbin_i<-elevation_lastbin[shot_number_id]
  z=rev(seq(elevation_lastbin_i,elevation_bin0_i,(elevation_bin0_i-elevation_lastbin_i)/rx_sample_count[shot_number_id]))[-1]

  waveform<-new("gedi.fullwaveform", dt = data.table::data.table(rxwaveform=rxwaveform_i,elevation=z))

  return(waveform)
 }
}


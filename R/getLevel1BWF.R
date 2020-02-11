#'Get GEDI Level1B waveform
#'
#'@description Extract GEDI Level1B waveform for a given shot number
#'
#'
#'@param level1b h5file; S4 object of class H5File
#'@param shot_number Shot number, an vector represeting
#'
#'@return Returns An object of class "gedi.waveform";
#'
#'@examples
#'
#'#' LVIS level 2 file path
#'level1_filepath = system.file("extdata", "lvis_level1_clip.h5", package="rLVIS")
#'
#'# Rectangle
#'xleft = 9.35986
#'xright = 9.35988
#'ybottom = 0.5786
#'ytop = 0.5790
#'
#'#' Reading LVIS level 2 file
#'level1_waveform = readLevel1b(level1_filepath)
#'
#'output = tempfile(fileext="h5")
#'
#'clipped_waveform = clipLevel1(level1_waveform, output, xleft, xright, ybottom, ytop)
#'
#'@export
#'
getxWF<-function(x,shot_number){
  x<-x@h5
  groups_id<-grep("BEAM\\d{4}$",gsub("/","",
                                     list.groups(x, recursive = F)), value = T)

  i = NULL
  #k<-"BEAM1011"
  for ( k in groups_id){
    gid<-max(x[[paste0(k,"/shot_number")]][]==shot_number)
    if (gid==1) {i=k}
  }

  if(is.null(i)) {
    stop(paste0("Shot number ", shot_number, " was not found within the dataset!"))
  }

  shot_number_i<-x[[paste0(i,"/shot_number")]][]
  shot_number_id<-which(shot_number_i[]==shot_number)
  elevation_bin0<-x[[paste0(i,"/geolocation/elevation_bin0")]][]
  elevation_lastbin<-x[[paste0(i,"/geolocation/elevation_lastbin")]][]
  rx_sample_count<-x[[paste0(i,"/rx_sample_count")]][]
  rx_sample_start_index<-x[[paste0(i,"/rx_sample_start_index")]][]
  rxwaveform_i<-x[[paste0(i,"/rxwaveform")]][rx_sample_start_index[shot_number_id]:(rx_sample_start_index[shot_number_id]+rx_sample_count[shot_number_id]-1)]
  rxwaveform_inorm<-(rxwaveform_i-min(rxwaveform_i))/(max(rxwaveform_i)-min(rxwaveform_i))*100
  elevation_bin0_i<-elevation_bin0[shot_number_id]
  elevation_lastbin_i<-elevation_lastbin[shot_number_id]
  z=rev(seq(elevation_lastbin_i,elevation_bin0_i,(elevation_bin0_i-elevation_lastbin_i)/rx_sample_count[shot_number_id]))[-1]

  waveform<-new("gedi.fullwaveform", dt = data.table::data.table(rxwaveform=rxwaveform_i,elevation=z))

  return(waveform)
}


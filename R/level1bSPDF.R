#'Read LVIS Level1b data
#'
#'@description This function reads LVIS level1b data
#'
#'@usage readLevel1b(level1bpath)
#'
#'@param level1bpath file path pointing to LVIS level1b data (H5 format)
#'@return H5File; S4 object of class H5File;
#'@author Carlos Alberto Silva. This function calls \emph{h5file} function from h5 package (Author: Mario Annau)
#'@seealso \code{\link[h5]{h5file}} in the \emph{h5} package.
#'@examples
#'
#'# LVIS level1b file path
#'level1b_filepath_zip <- system.file("extdata", "LVIS_Mondah_level1b.zip", package="rLVIS")
#'unzip(level1b_filepath_zip, exdir = tempdir())
#'level1b_filepath <- file.path(tempdir(), "LVIS_Mondah_level1b.h5")
#'
#'# Reading LVIS level1b file
#'level1b<-readLevel1b(level1bpath=level1b_filepath)
#'
#'
#'@importFrom hdf5r H5File
#'@importFrom data.table data.table
#'@export

level1bSPDF<-function(level1b) {
    groups_id<-grep("BEAM\\d{4}$",gsub("/","",
                                     list.groups(level1b, recursive = F)), value = T)
    pb <- txtProgressBar(min = 0, max = length(groups_id), style = 3)
    level1b.dt<-data.table()
    i.s=0
    #i="BEAM0010"
    for ( i in gsub("/","",groups_id)){
        i.s<-i.s+1
        setTxtProgressBar(pb, i.s)
        level1b.dt<-rbind(level1b.dt, cbind(
        latitude_bin0=level1b[[paste0(i,"/geolocation/latitude_bin0")]][],
        latitude_lastbin=level1b[[paste0(i,"/geolocation/latitude_lastbin")]][],
        longitude_bin0=level1b[[paste0(i,"/geolocation/longitude_bin0")]][],
        longitude_lastbin=level1b[[paste0(i,"/geolocation/longitude_lastbin")]][],
        elevation_bin0=level1b[[paste0(i,"/geolocation/elevation_bin0")]][],
        elevation_lastbin=level1b[[paste0(i,"/geolocation/elevation_lastbin")]][],
        beam=level1b[[paste0(i,"/beam")]][],
        shot_number=as.character(level1b[[paste0(i,"/shot_number")]][])))
  }
  level1b.spdf<-SpatialPointsDataFrame(cbind(as.numeric(level1b.dt$latitude_bin0),as.numeric(level1b.dt$longitude_bin0)),data=level1b.dt)
  level1b.spdf<- new("gedi.level1bSPDF", spdf = level1b.spdf)
  close(pb)
  return(level1b.spdf)
}

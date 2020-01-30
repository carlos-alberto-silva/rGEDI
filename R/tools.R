<<<<<<< HEAD
#'Read LVIS Level1b data
#'
#'@description This function reads LVIS level1b data
#'
#'
#'@param level1b H5File object from level1b
#'@return H5File; S4 object of class H5File;
#'@author Carlos Alberto Silva. This function calls \emph{h5file} function from h5 package (Author: Mario Annau)
#'@seealso \code{\link[h5]{h5file}} in the \emph{h5} package.
#'@importFrom hdf5r H5File
#'@importFrom data.table data.table
#'@export
level1b.dt<-function(level1b) {
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
  close(pb)
  return(level1b.dt)
}
=======
#'Read LVIS Level1b data
#'
#'@description This function reads LVIS level1b data
#'
#'
#'@param level1b H5File object from level1b
#'@return H5File; S4 object of class H5File;
#'@author Carlos Alberto Silva. This function calls \emph{h5file} function from h5 package (Author: Mario Annau)
#'@seealso \code{\link[h5]{h5file}} in the \emph{h5} package.
#'@import bit64
#'@importFrom hdf5r H5File
#'@importFrom data.table data.table
#'@export
level1b.dt<-function(level1b) {
    groups_id<-grep("BEAM\\d{4}$",gsub("/","",
                                     list.groups(level1b, recursive = F)), value = T)
    pb <- txtProgressBar(min = 0, max = length(groups_id), style = 3)
    level1b.dt<-data.table()
    i.s=0
    #i="BEAM0010"
    for ( i in gsub("/","",groups_id)){
        i.s<-i.s+1
        setTxtProgressBar(pb, i.s)
        shot_number.df <- bit64::as.data.frame.integer64(level1b[[paste0(i,"/shot_number")]][])
        names(shot_number.df) <- "shot_number"
        level1b.dt <- rbind(level1b.dt, cbind(
          latitude_bin0=level1b[[paste0(i,"/geolocation/latitude_bin0")]][],
          latitude_lastbin=level1b[[paste0(i,"/geolocation/latitude_lastbin")]][],
          longitude_bin0=level1b[[paste0(i,"/geolocation/longitude_bin0")]][],
          longitude_lastbin=level1b[[paste0(i,"/geolocation/longitude_lastbin")]][],
          elevation_bin0=level1b[[paste0(i,"/geolocation/elevation_bin0")]][],
          elevation_lastbin=level1b[[paste0(i,"/geolocation/elevation_lastbin")]][],
          beam=level1b[[paste0(i,"/beam")]][],
          shot_number.df))
  }
  close(pb)
  return(level1b.dt)
}
>>>>>>> f57d0feba7a47ee64476a94c1d808eb822c65e80

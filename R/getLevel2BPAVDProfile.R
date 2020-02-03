#'Read LVIS level2a data
#'
#'@description This function reads LVIS level2a data
#'
#'
#'@param level2a H5File object from level2a
#'@return H5File; S4 object of class H5File;
#'@author Carlos Alberto Silva. This function calls \emph{h5file} function from h5 package (Author: Mario Annau)
#'@seealso \code{\link[h5]{h5file}} in the \emph{h5} package.
#'
#'@import hdf5r
#'@import utils
#'@import sp
#'@importFrom hdf5r H5File
#'@importFrom data.table data.table
#'@importFrom sp SpatialPointsDataFrame
#'@export
getLevel2BPAVDProfile<-function(level2b){
  level2b<-level2b@data
  groups_id<-grep("BEAM\\d{4}$",gsub("/","",
                                     hdf5r::list.groups(level2b, recursive = F)), value = T)
  m.dt<-data.table::data.table()
  pb <- utils::txtProgressBar(min = 0, max = length(groups_id), style = 3)
  i.s=0
  for ( i in groups_id){
    i.s<-i.s+1
    utils::setTxtProgressBar(pb, i.s)
    level2b_i<-level2b[[i]]
    m<-data.table::data.table(
      beam<-rep(i,length(level2b_i[["shot_number"]][])),
      shot_number=level2b_i[["shot_number"]][],
      lat_lowestmode=level2b_i[["geolocation/lat_lowestmode"]][],
      lon_lowestmode=level2b_i[["geolocation/lon_lowestmode"]][],
      elev_highestreturn=level2b_i[["geolocation/elev_highestreturn"]][],
      elev_lowestmode=level2b_i[["geolocation/elev_lowestmode"]][],
      height_lastbin=level2b_i[["geolocation/height_lastbin"]][],
      height_bin0=level2b_i[["geolocation/height_bin0"]][],
      pavd_z=t(level2b_i[["pavd_z"]][,1:level2b_i[["pavd_z"]]$dims[2]]))
    m.dt<-rbind(m.dt,m)
  }
  colnames(m.dt)<-c("beam","shot_number","lat_lowestmode",
                    "lon_lowestmode","elev_highestreturn",
                    "elev_lowestmode","height_lastbin",
                    "height_bin0",paste0("pavd_z",seq(5,30*5,5),"m"))
  close(pb)
  return(m.dt)
}

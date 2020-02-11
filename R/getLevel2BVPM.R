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
#'
#'
getxVPM<-function(x){
  x<-x@h5
  groups_id<-grep("BEAM\\d{4}$",gsub("/","",
                                     hdf5r::list.groups(x, recursive = F)), value = T)
  m.dt<-data.table::data.table()
  pb <- utils::txtProgressBar(min = 0, max = length(groups_id), style = 3)
  i.s=0
  for ( i in groups_id){
    i.s<-i.s+1
    utils::setTxtProgressBar(pb, i.s)
    x_i<-x[[i]]
    m<-data.table::data.table(
      beam<-rep(i,length(x_i[["shot_number"]][])),
      shot_number=x_i[["shot_number"]][],
      delta_time=x_i[["delta_time"]][],
      lat_lowestmode=x_i[["geolocation/lat_lowestmode"]][],
      lon_lowestmode=x_i[["geolocation/lon_lowestmode"]][],
      elev_highestreturn=x_i[["geolocation/elev_highestreturn"]][],
      elev_lowestmode=x_i[["geolocation/elev_lowestmode"]][],
      pai=x_i[["pai"]][],
      fhd_normal=x_i[["fhd_normal"]][],
      omega=x_i[["omega"]][],
      pgap_theta=x_i[["pgap_theta"]][],
      cover=x_i[["cover"]][])
    m.dt<-rbind(m.dt,m)
  }
  colnames(m.dt)<-c("beam","shot_number","delta_time","lat_lowestmode","lon_lowestmode",
                    "elev_highestreturn","elev_lowestmode","pai",
                    "fhd_normal","omega","pgap_theta","cover")
  close(pb)
  return(m.dt)
}



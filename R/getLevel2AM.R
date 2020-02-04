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
getLevel2AM<-function(level2a){
  level2a<-level2a@h5
  groups_id<-grep("BEAM\\d{4}$",gsub("/","",
                                     hdf5r::list.groups(level2a, recursive = F)), value = T)
  rh.dt<-data.table::data.table()
  pb <- utils::txtProgressBar(min = 0, max = length(groups_id), style = 3)
  i.s=0
  for ( i in groups_id){
    i.s<-i.s+1
    utils::setTxtProgressBar(pb, i.s)
    level2a_i<-level2a[[i]]
    rhs<-data.table::data.table(
      beam<-rep(i,length(level2a_i[["shot_number"]][])),
      shot_number=level2a_i[["shot_number"]][],
      lat_lowestmode=level2a_i[["lat_lowestmode"]][],
      lon_lowestmode=level2a_i[["lon_lowestmode"]][],
      elev_highestreturn=level2a_i[["elev_highestreturn"]][],
      elev_lowestmode=level2a_i[["elev_lowestmode"]][],
      t(level2a_i[["rh"]][,level2a_i[["rh"]]$dims[2]]))
    rh.dt<-rbind(rh.dt,rhs)
    }
  colnames(rh.dt)<-c("beam","shot_number","lat_lowestmode","lon_lowestmode",
                     "elev_highestreturn","elev_lowestmode",paste0("rh",seq(0,100)))
  close(pb)
  return(rh.dt)
}



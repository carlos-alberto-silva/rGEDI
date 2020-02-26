#'Get GEDI Plant Area Volume Density (PAVD) Index (GEDI Level2B)
#'
#'@description This function extracts Plant Area Volume Density (PAVD) Index from GEDI Level2B data.
#'
#'@usage getLevel2BPAVDProfile(level2b)
#'
#'@param level2b A GEDI Level2B object (output of \code{\link[rGEDI:readLevel2B]{readLevel2B}} function). A S4 object of class "gedi.level2b".
#'
#'@return An S4 object of class \code{\link[data.table:data.table]{data.table-class}}
#'containing the PAVD Index.
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi02_bv001/
#'
#'@examples
#'# specify the path to GEDI Level 2B data
#'level2bpath <- system.file("extdata", "GEDIexample_level01B.h5", package="rGEDI")
#'
#'# Reading GEDI level2B data
#'level2b <- readLevel2B(level2bpath)
#'
#'# Getting GEDI Plant Area Volume Density (PAVD) Index
#'level2BPAVDProfile<-getLevel2BPAVDProfile(level2b)
#'head(level2BPAVDProfile)
#'
#'@export
getLevel2BPAVDProfile<-function(level2b){
  level2b<-level2b@h5
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
                    "height_bin0",paste0("pavd_z",seq(0,30*5,5)[-31],"_",seq(5,30*5,5),"m"))
  close(pb)
  return(m.dt)
}

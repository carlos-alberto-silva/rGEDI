#'Get GEDI Elevation and Height Metrics (GEDI Level2A)
#'
#'@description This function extracts Elevation and Relative Height (RH) metrics from GEDI Level2A data .
#'
#'@usage getLevel2AM(level2a)
#'
#'@param level2a A GEDI Level2A object (output of \code{\link[rGEDI:readLevel2A]{readLevel2A}} function). A S4 object of class "gedi.level2a".
#'
#'@return An S4 object of class \code{\link[data.table:data.table]{data.table-class}}
#'containing the elevation and relative heights.
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi02_av001/
#'
#'@examples
#'
#'# specify the path to GEDI level2A data (zip file)
#'level2A_fp_zip <- system.file("extdata",
#'                   "GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level2A data
#'level2Apath <- unzip(level2A_fp_zip,exdir = dirname(level2A_fp_zip))
#'
#'# Reading GEDI level2A data (h5 file)
#'level2a<-readLevel2A(level2Apath=level2Apath)
#'
#'# Get GEDI Elevation and Height Metrics
#'level2AM<-getLevel2AM(level2a)
#'head(level2AM)
#'
#'close(level2a)
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

    if (any(hdf5r::list.datasets(level2a_i)=="shot_number")){

    if(length(level2a_i[["rh"]]$dims)==2) {
      rh=t(level2a_i[["rh"]][,])
    } else {
      rh=t(level2a_i[["rh"]][])
    }

    rhs<-data.table::data.table(
      beam<-rep(i,length(level2a_i[["shot_number"]][])),
      shot_number=level2a_i[["shot_number"]][],
      lat_lowestmode=level2a_i[["lat_lowestmode"]][],
      lon_lowestmode=level2a_i[["lon_lowestmode"]][],
      elev_highestreturn=level2a_i[["elev_highestreturn"]][],
      elev_lowestmode=level2a_i[["elev_lowestmode"]][],
      rh)
    rh.dt<-rbind(rh.dt,rhs)
    }
  }

  colnames(rh.dt)<-c("beam","shot_number","lat_lowestmode","lon_lowestmode",
                     "elev_highestreturn","elev_lowestmode",paste0("rh",seq(0,100)))
  close(pb)
  return(rh.dt)
}



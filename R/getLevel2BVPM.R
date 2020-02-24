#'Get Vegetation Profile Biophysical Variables (GEDI Level2B)
#'
#'@description This function extracts information from GEDI Level2B data:
#'Total Plant Area Index,	Foliage Height Diversity, Foliage Clumping Index,
#'Total Gap Probability (theta), and Total canopy cover.
#'
#'@usage getLevel2BVPM(level2b)
#'
#'@param level2b A GEDI Level2B object (output of \code{\link[rGEDI:readLevel2B]{readLevel2B}}
#'function). A S4 object of class "gedi.level2b".
#'
#'#'@return An S4 object of class \code{\link[data.table:data.table]{data.table-class}}
#'containing the Vegetation Profile Biophysical Variables.
#'
#'@details These are the biophysical variables extracted:
#'\itemize{
#'\item \emph{pai} - Total Plant Area Index.
#'\item \emph{fhd_normal} -	Foliage Height Diversity.
#'\item \emph{omega} -	Foliage Clumping Index.
#'\item \emph{pgap_theta} -	Total Gap Probability (theta).
#'\item \emph{cover} -	Total canopy cover.
#'}
#'
#'@examples
#'# specify the path to GEDI Level 2B data
#'level2bpath <- system.file("extdata", "GEDIexample_level01B.h5", package="rGEDI")
#'
#'# Reading GEDI level2B data
#'level2b <- readLevel2B(level2bpath)
#'
#'# Get GEDI level2B geolocations
#'level2BVPM<-getLevel2BVPM(level2b)
#'head(level2BVPM)
#'
#'@export
getLevel2BVPM<-function(level2b){
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
      delta_time=level2b_i[["delta_time"]][],
      lat_lowestmode=level2b_i[["geolocation/lat_lowestmode"]][],
      lon_lowestmode=level2b_i[["geolocation/lon_lowestmode"]][],
      elev_highestreturn=level2b_i[["geolocation/elev_highestreturn"]][],
      elev_lowestmode=level2b_i[["geolocation/elev_lowestmode"]][],
      pai=level2b_i[["pai"]][],
      fhd_normal=level2b_i[["fhd_normal"]][],
      omega=level2b_i[["omega"]][],
      pgap_theta=level2b_i[["pgap_theta"]][],
      cover=level2b_i[["cover"]][])
    m.dt<-rbind(m.dt,m)
  }
  colnames(m.dt)<-c("beam","shot_number","delta_time","lat_lowestmode","lon_lowestmode",
                    "elev_highestreturn","elev_lowestmode","pai",
                    "fhd_normal","omega","pgap_theta","cover")
  close(pb)
  return(m.dt)
}



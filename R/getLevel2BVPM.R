#'Get GEDI Vegetation Profile Biophysical Variables (GEDI Level2B)
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
#'@seealso https://lpdaac.usgs.gov/products/gedi02_bv001/
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
#'\donttest{
#'# specify the path to download GEDI example dataset
#'outdir<-getwd()
#'
#'# downloading GEDI example dataset (zip file)
#'download.file(
#'              paste0(
#'                     "https://github.com/carlos-alberto-silva/rGEDI/",
#'                     "releases/download/examples/examples.zip"
#'              ),
#'              destfile=file.path(outdir, "examples.zip"))
#'
#'# unzip the file
#'unzip(file.path(outdir, "examples.zip"))
#'
#'# specify the path to GEDI level2B data
#'level2bpath = file.path(outdir, "GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.h5")
#'
#'# Reading GEDI level1B file
#'level2b<-readLevel2b(gedilevel2b)
#'
#'# Get GEDI Vegetation Profile Biophysical Variables
#'level2BVPM<-getLevel2BVPM(level2b)
#'head(level2BVPM)
#'
#'close(level2b)
#'}
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
      latitude_lastbin=level2b_i[["geolocation/latitude_lastbin"]][],
      latitude_bin0=level2b_i[["geolocation/latitude_bin0"]][],
      longitude_bin0=level2b_i[["geolocation/longitude_bin0"]][],
      longitude_lastbin=level2b_i[["geolocation/longitude_lastbin"]][],
      elev_highestreturn=level2b_i[["geolocation/elev_highestreturn"]][],
      elev_lowestmode=level2b_i[["geolocation/elev_lowestmode"]][],
      pai=level2b_i[["pai"]][],
      fhd_normal=level2b_i[["fhd_normal"]][],
      omega=level2b_i[["omega"]][],
      pgap_theta=level2b_i[["pgap_theta"]][],
      cover=level2b_i[["cover"]][])
    m.dt<-rbind(m.dt,m)
  }
  colnames(m.dt)<-c("beam","shot_number","delta_time","latitude_lastbin","latitude_bin0",
                    "longitude_lastbin","longitude_bin0",
                    "elev_highestreturn","elev_lowestmode","pai",
                    "fhd_normal","omega","pgap_theta","cover")
  close(pb)
  return(m.dt)
}

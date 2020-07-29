
#'Get GEDI Canopy Cover and Vertical Profile Metrics (GEDI Level2B)
#'
#'@description This function extracts information from GEDI Level2B data:
#'Total Plant Area Index,	Foliage Height Diversity, Foliage Clumping Index,
#'Total Gap Probability (theta), and Total canopy cover.
#'
#'@usage getLevel2BVPM(level2b)
#'
#'@param level2b A GEDI Level2B object (output of \code{\link[rGEDI:readLevel2B]{readLevel2B}}
#'function). An S4 object of class "gedi.level2b".
#'
#'@return Returns an S4 object of class \code{\link[data.table:data.table]{data.table-class}}
#'containing the Vegetation Profile Biophysical Variables.
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi02_bv001/
#'
#'@details These are the biophysical variables and additional information extracted:
#'\itemize{
#'\item \emph{beam} Beam identifie
#'\item \emph{shot_number} Shot number
#'\item \emph{algorithmrun_flag} The L2B algorithm is run if this flag is set to 1 indicating data have sufficient waveform fidelity for L2B to run
#'\item \emph{l2b_quality_flag} L2B quality flag
#'\item \emph{delta_time} Transmit time of the shot since Jan 1 00:00 2018
#'\item \emph{sensitivity} Maxmimum canopy cover that can be penetrated
#'\item \emph{solar_elevation} Solar elevation
#'\item \emph{latitude_lastbin} Latitude of last bin of the pgap_theta_z, interpolated from L1B waveform coordinate
#'\item \emph{latitude_bin0} Latitude of first bin of the pgap_theta_z, interpolated from L1B waveform coordinate
#'\item \emph{elev_highestreturn} Elevation of highest detected return relative to reference ellipsoid
#'\item \emph{elev_lowestmode} Elevation of center of lowest mode relative to reference ellipsoid
#'\item \emph{pai} Total Plant Area Index
#'\item \emph{fhd_normal} Foliage Height Diversity
#'\item \emph{omega} Foliage Clumping Index
#'\item \emph{pgap_theta} Total Gap Probability (theta)
#'\item \emph{cover} Total canopy cover
#'}
#'
#'@examples
#'# Specifying the path to GEDI level2B data (zip file)
#'outdir = tempdir()
#'level2B_fp_zip <- system.file("extdata",
#'                   "GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level2A data
#'level2Bpath <- unzip(level2B_fp_zip,exdir = outdir)
#'
#'# Reading GEDI level2B data (h5 file)
#'level2b<-readLevel2B(level2Bpath=level2Bpath)
#'
#'# Extracting GEDI Vegetation Profile Biophysical Variables
#'level2BVPM<-getLevel2BVPM(level2b)
#'head(level2BVPM)
#'
#'close(level2b)
#'@export
getLevel2BVPM<-function(level2b){
  level2b<-level2b@h5
  groups_id<-grep("BEAM\\d{4}$",gsub("/","",
                                     hdf5r::list.groups(level2b, recursive = F)), value = T)
  m.dt<-data.table::data.table()
  pb <- utils::txtProgressBar(min = 0, max = length(groups_id), style = 3)
  i.s=0
  var.map = data.table::data.table(t(data.frame(list(
    # COL_NAMES              # H5_ADDRESS
    c("shot_number",         "shot_number"),
    c("algorithmrun_flag",   "algorithmrun_flag"),
    c("l2b_quality_flag",    "l2b_quality_flag"),
    c("delta_time",          "geolocation/delta_time"),
    c("sensitivity",         "sensitivity"),
    c("solar_elevation",     "geolocation/solar_elevation"),
    c("latitude_lastbin",    "geolocation/latitude_lastbin"),
    c("latitude_bin0",       "geolocation/latitude_bin0"),
    c("longitude_bin0",      "geolocation/longitude_bin0"),
    c("longitude_lastbin",   "geolocation/longitude_lastbin"),
    c("elev_highestreturn",  "geolocation/elev_highestreturn"),
    c("elev_lowestmode",     "geolocation/elev_lowestmode"),
    c("rh100",               "rh100"),
    c("pai",                 "pai"),
    c("fhd_normal",          "fhd_normal"),
    c("omega",               "omega"),
    c("pgap_theta",          "pgap_theta"),
    c("cover",               "cover")
  ))))
  colnames(var.map) = c("COL_NAMES", "H5_ADDRESS")

  for ( i in groups_id){
    i.s<-i.s+1
    utils::setTxtProgressBar(pb, i.s)
    level2b_i<-level2b[[i]]
    m<-data.table::data.table(
      beam=rep(i,length(level2b_i[["shot_number"]][])))
    for (row_index in 1:nrow(var.map)) {
      colname = var.map$COL_NAMES[row_index]
      h5.address = var.map$H5_ADDRESS[row_index]
      m[[colname]] <- level2b_i[[h5.address]][]
    }
    m.dt<-rbind(m.dt,m)
  }
  colnames(m.dt)<-c("beam", var.map$COL_NAMES)
  close(pb)
  return(m.dt)
}

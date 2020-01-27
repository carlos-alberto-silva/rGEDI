#'Read LVIS Level2 data
#'
#'@description This function reads LVIS level2 data
#'
#'@usage readLevel2(level2path, spdf=TRUE, glatlon=TRUE)
#'
#'@param level2path file path pointing to LVIS level2 data (txt format)
#'@param spdf if true, outups an object of class \code{SpatialPointsDataFrame}
#'@param glatlon if true, GLON and GLAT will be used for creating the \code{SpatialPointsDataFrame}.
#'If false, TLON and TLAT will be used instead
#'@return An object of class \code{SpatialPointsDataFrame} or \code{data.table};
#'@author Carlos Alberto Silva.
#'@examples
#'
#'# LVIS level2 file path
#'level2_filepath_zip <- system.file("extdata", "LVIS_Mondah_level2.zip", package="rLVIS")
#'unzip(level2_filepath_zip, exdir = tempdir())
#'level2_filepath <- file.path(tempdir(), "LVIS_Mondah_level2.txt")
#'
#'# Reading LVIS level1b file
#'level2_spdf<-readLevel2(level2path=level2_filepath, spdf=TRUE, glatlon=TRUE)
#'
#'#' Plot LVIS Level2 data
#'plotLevel2(level2_spdf=level2_spdf, color = "RH100",
#'            colorPalette = c("blue","green","yellow","red"))
#'
#'
#'@export
readLevel2<-function(level2path, spdf=TRUE, glatlon=TRUE) {

  l2file<-utils::read.table(level2path,sep="")
  colnames(l2file)<-c("LFID","SHOTNUMBER","TIME","GLON","GLAT","ZG","TLON","TLAT","ZT","RH10","RH15","RH20","RH25",
                      "RH30","RH35","RH40","RH45","RH50","RH55","RH60","RH65","RH70","RH75","RH80","RH85",
                      "RH90","RH95","RH96","RH97","RH98","RH99","RH100","AZIMUTH","INCIDENTANGLE","RANGE","FLAG1",
                      "FLAG2","FLAG3")
  if (spdf==TRUE){

    if (glatlon==TRUE) {
      Level2<-sp::SpatialPointsDataFrame(l2file[,c("GLON", "GLAT")],data=l2file)
      sp::proj4string(Level2) <- sp::CRS("+proj=longlat +datum=WGS84")

    } else {

      Level2<-sp::SpatialPointsDataFrame(l2file[,c("TLON", "TLAT")],data=l2file)
      sp::proj4string(Level2) <- sp::CRS("+proj=longlat +datum=WGS84")

    }
  } else {

    Level2<-data.table::data.table(l2file)

  }
  return(Level2)
}

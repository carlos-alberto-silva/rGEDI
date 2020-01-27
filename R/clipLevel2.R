#'Clip LVIS Level2 data
#'
#'@description Clip LVIS Level2 data within a given geometry
#'
#'@usage clipLevel2(level2_spdf, polygon_spdf)
#'
#'@param level2_spdf h5file; S4 object of class H5File
#'@param polygon_spdf dataframe containing LVIS level 2 data
#'@return Returns An object of class \code{SpatialPoligonDataFrame} ; subset of LVIS Level2 data
#'@author Carlos Alberto Silva.
#'@examples
#'
#'#' LVIS level 2 file path
#'level2_filepath_zip <- system.file("extdata", "LVIS_Mondah_level2.zip", package="rLVIS")
#'unzip(level2_filepath_zip, exdir = tempdir())
#'level2_filepath <- file.path(tempdir(), "LVIS_Mondah_level2.txt")
#'
#'#' Polgons file path
#'polygons_filepath <- system.file("extdata", "LVIS_Mondah_polygons.shp", package="rLVIS")
#'
#'#' Reading LVIS level 2 file
#'level2_spdf<-readLevel2(level2path=level2_filepath,spdf=TRUE,glatlon=TRUE)
#'
#'#' Plot LVIS Level2 data
#'plotLevel2(level2_spdf=level2_spdf, color = "RH100",
#'           colorPalette = c("blue","green","yellow","red"))
#'
#'#' Reading Polygons
#'library(rgdal)
#'plots<-readOGR(polygons_filepath)
#'proj4string(plots) <- CRS("+proj=longlat +datum=WGS84")
#'plot(plots, add=TRUE, border="black", lwd=2)
#'
#'#'Clipping LVIS Level2 data
#'level2_spdf_sub<-clipLevel2(level2_spdf=level2_spdf,polygon_spdf=plots)
#'
#'#' Plot LVIS Level2 data
#'plotLevel2(level2_spdf=level2_spdf_sub, color = "RH100",
#'          colorPalette = c("blue","green","yellow","red"))
#'
#'@export
clipLevel2<-function(level2_spdf, polygon_spdf){
  CLIPID<-sp::over(level2_spdf,polygon_spdf)
  level2_spdf@data<-cbind(level2_spdf@data,CLIPID=CLIPID[,1])
  level2_spdf<-level2_spdf[rownames(CLIPID)[!is.na(CLIPID)],]
  return(level2_spdf)
}

#'LVIS Level2 data visualization in 2-D
#'
#'@description This function plots LVIS Level2 in 2-D
#'
#'@usage plotLevel2(level2_spdf, color, colorPalette,...)
#'
#'@param level2_spdf LVIS l2 dataset; object of class \code{SpatialPointsDataFrame}
#'@param color  The field name used to color the points.  Default is RH100
#'@param colorPalette A vector defining the color  palette.
#'@param ... passing arguments on to the plot function
#'@return A 2-D figure of LVIS Level2 data;
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
#'head(level2_spdf@data)
#'plotLevel2(level2_spdf=level2_spdf, color = "RH100",
#'           colorPalette = c("blue","green","yellow","red"),
#'           axes=TRUE)
#'grid()
#'@export
plotLevel2 = function(level2_spdf, color = "RH100", colorPalette = c("blue","green","yellow","red"),...)
{
  v <- (level2_spdf@data[,color] - min(level2_spdf@data[,color]))/diff(range(level2_spdf@data[,color]))
  x <- colorRamp(colorPalette)(v)
  col<-rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
  sp::plot(level2_spdf, col=col,...)
}

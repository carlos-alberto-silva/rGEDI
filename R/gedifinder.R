#'GEDI filder
#'
#'@description This function finds the path to GEDI data within a boundary box provided
#'
#'
#'@param x GEDI data type; "GEDI01_B", "GEDI02_A" or "GEDI02_B"
#'@param xleft numeric. left x coordinates of rectangles (degree).
#'@param xright numeric. right x coordinates of rectangles.
#'@param ybottom numeric. bottom y coordinates of rectangles.
#'@param ytop numeric. top y coordinates of rectangles.
#'@return vector point out the path for download GEDI data within the bbox provided
#'@seealso gedifilder tool: https://lpdaacsvc.cr.usgs.gov/services/gedifinder
#'
#'@examples
#'
#'gedi02_blist<-GEDIfinder(x="GEDI02_B",xleft=-40.9218,xright=-42.0248,ybottom=-71.7736,ytop=-73.7731)
#'
#'@export
GEDIfinder<-function(x, xleft, xright, ybottom, ytop){

  response = httr::GET(paste0("https://lpdaacsvc.cr.usgs.gov/services/gedifinder?product=",x,"&version=001&bbox=",xleft,",",ytop,",",xright,",",ybottom,"&output=html"))
  content = httr::content(response, "text")
  links = gsub(".*a href=([^<]*).*", "\\1", strsplit(content,"<br/>")[[1]])
  return(links)

}



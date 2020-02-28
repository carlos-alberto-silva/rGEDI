#'GEDI finder
#'
#'@description This function finds the path to GEDI data within a boundary box coordinates provided
#'
#'#'@usage gedifinder(level2BPAVDProfile, xleft, xright, ybottom, ytop)
#'
#'@param level GEDI data level; Options: "GEDI01_B", "GEDI02_A" or "GEDI02_B"
#'@param xmin Numeric. West longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param xmax Numeric. East longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param ymin Numeric. South latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'@param ymax Numeric. North latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'
#'@return vector object pointing out the path for download GEDI data within
#'the boundary box coordinates provided
#'@seealso This function just call the gedifilder tool developted by LPDAAC: https://lpdaacsvc.cr.usgs.gov/services/gedifinder
#'
#'@examples
#'
#'# specify bounding box coordinates
#'xmin<- -44.17246
#'ymin<- -44.0654
#'xmax<- -13.76913
#'ymax<- -13.67646
#'
#'# Getting the path to GEDI data for the specified boundary box coordinates
#'gedi02b_list<-gedifinder(level="GEDI02_B",xleft=xleft,xright=xright,ybottom=ybottom,ytop=ytop)
#'
#'@export
gedifinder<-function(level, xmin, xmax, ymin, ymax){
  response = httr::GET(paste0("https://lpdaacsvc.cr.usgs.gov/services/gedifinder?product=",level,"&version=001&bbox=",xmin,",",ymax,",",xmax,",",ymin,"&output=html"))
  content = httr::content(response, "text")
  links = gsub(".*a href=([^<]*).*", "\\1", strsplit(content,"<br/>")[[1]])
  return(links)
}



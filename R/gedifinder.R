#'GEDI finder
#'
#'@description This function finds the path to GEDI data within a boundary box coordinates provided
#'
#'#'@usage gedifinder(level2BPAVDProfile, xleft, xright, ybottom, ytop, output="")
#'
#'@param level GEDI data level; Options: "GEDI01_B", "GEDI02_A" or "GEDI02_B"
#'@param xleft numeric. left x coordinates of rectangles (degree).
#'@param xright numeric. right x coordinates of rectangles (degree).
#'@param ybottom numeric. bottom y coordinates of rectangles (degree).
#'@param ytop numeric. top y coordinates of rectangles (degree).
#'@return vector object pointing out the path for download GEDI data within
#'the boundary box coordinates provided
#'@seealso This function just call the gedifilder tool developted by LPDAAC: https://lpdaacsvc.cr.usgs.gov/services/gedifinder
#'
#'@examples
#'
#'# Bounding rectangle coordinates
#'xleft = -40.9218
#'xright = -42.0248
#'ybottom = -71.7736
#'ytop = -73.7731
#'
#'# Getting the path to GEDI data within the provided boundary box coordinates
#'gedi02_b_list<-gedifinder(level="GEDI02_B",xleft=xleft,xright=xright,ybottom=ybottom,ytop=ytop)
#'
#'@export
gedifinder<-function(level, xleft, xright, ybottom, ytop){
  response = httr::GET(paste0("https://lpdaacsvc.cr.usgs.gov/services/gedifinder?product=",level,"&version=001&bbox=",xleft,",",ytop,",",xright,",",ybottom,"&output=html"))
  content = httr::content(response, "text")
  links = gsub(".*a href=([^<]*).*", "\\1", strsplit(content,"<br/>")[[1]])
  return(links)
}



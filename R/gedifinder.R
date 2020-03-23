#'GEDI finder
#'
#'@description This function finds the path to GEDI data within a boundary box coordinates provided
#'
#'#'@usage gedifinder(level="GEDI02_B",xmin,xmax,ymin,ymax)
#'
#'@param level GEDI data level; Options: "GEDI01_B", "GEDI02_A" or "GEDI02_B"
#'@param ul_lat Numeric. Upper left (ul) corner coordinates, in lat (decimal degrees) for the bounding box of the area of interest.
#'@param ul_lon Numeric. Upper left (ul) corner coordinates, in lon (decimal degrees) for the bounding box of the area of interest.
#'@param lr_lat Numeric. Lower right (ul) corner coordinates, in lat (decimal degrees) for the bounding box of the area of interest.
#'@param lr_lon Numeric. Lower right (ul) corner coordinates, in lon (decimal degrees) for the bounding box of the area of interest.
#'
#'@return Return a vector object pointing out the path saving the downloaded GEDI data within
#'the boundary box coordinates provided
#'@seealso bbox: The correct format is upper left and lower right corner coordinates, in lat,lon ordering, for the bounding box of the area of interest (e.g. [ul_lat,ul_lon,lr_lat,lr_lon]).
#'This function just call the gedifinder tool developted by LPDAAC:
#'https://lpdaacsvc.cr.usgs.gov/services/gedifinder
#'
#'@examples
#'\donttest{
#' # gedifinder is a web service provided by NASA
#' # usually the request takes more than 5 seconds
#'
#'# Specifying bounding box coordinates
#'ul_lat<- 42.0
#'ul_lon<- -100
#'lr_lat<- 40.0
#'lr_lon<- -96.0
#'
#'# Extracting the path to GEDI data for the specified boundary box coordinates
#'gedi02b_list<-gedifinder(level="GEDI02_B",ul_lat, ul_lon, lr_lat, lr_lon)
#'}
#'@import jsonlite curl
#'@export
gedifinder<-function(level, ul_lat, ul_lon, lr_lat, lr_lon){
  response = curl::curl(sprintf(
    "https://lpdaacsvc.cr.usgs.gov/services/gedifinder?%s=%s&version=001&%s=%f,%f,%f,%f&output=json",
    "product", level,
    "bbox", ul_lat, ul_lon, lr_lat,lr_lon))
  content = suppressWarnings(readLines(response))
  close(response)
  results = jsonlite::parse_json(content)
  if(results$message != "") {
    stop(results$message)
  }
  return(simplify2array(results$data))
}

#'GEDI finder
#'
#'@description This function finds the exact granule(s) that contain GEDI data for a given region of interest and date range
#''
#'@param product GEDI data level; Options: "GEDI01_B", "GEDI02_A" or "GEDI02_B"
#'@param ul_lat Numeric. Upper left (ul) corner coordinates, in lat (decimal degrees) for the bounding box of the area of interest.
#'@param ul_lon Numeric. Upper left (ul) corner coordinates, in lon (decimal degrees) for the bounding box of the area of interest.
#'@param lr_lat Numeric. Lower right (ul) corner coordinates, in lat (decimal degrees) for the bounding box of the area of interest.
#'@param lr_lon Numeric. Lower right (ul) corner coordinates, in lon (decimal degrees) for the bounding box of the area of interest.
#'@param version Character. The version of the GEDI product files to be returned. Default "001".
#'@param daterange Vector. Date range. Specify your start and end dates (year-month-day). Ex.: c("2019-07-01","2020-05-22"). If NULL (default),
#'the date range filter will be not applied.
#'
#'@return Return a vector object pointing out the path saving the downloaded GEDI data within
#'the boundary box coordinates provided
#'@seealso bbox: Defined by the upper left and lower right corner coordinates, in lat,lon ordering, for the bounding box of the area of interest (e.g. [ul_lat,ul_lon,lr_lat,lr_lon]).
#'This function relies on the existing LP DAAC gedifinder tool:
#'https://lpdaacsvc.cr.usgs.gov/services/gedifinder
#'
#'@examples
#'\donttest{
#'# gedifinder is a web service provided by NASA
#'# usually the request takes more than 5 seconds
#'
#'# Specifying bounding box coordinates
#'ul_lat<- 42.0
#'ul_lon<- -100
#'lr_lat<- 40.0
#'lr_lon<- -96.0
#'
#'# Specifying the date range
#'daterange=c("2019-07-01","2020-05-22")
#'
#'# Extracting the path to GEDI data for the specified boundary box coordinates
#'gedi02b_list<-gedifinder(product="GEDI02_B",
#'                                 ul_lat,
#'                                 ul_lon,
#'                                 lr_lat,
#'                                 lr_lon,
#'                                 version="001",
#'                                 daterange=daterange)
#'
#'}
#'@import jsonlite curl
#'@export
gedifinder<-function(product, ul_lat, ul_lon, lr_lat, lr_lon,version="001",daterange=NULL){
  response = curl::curl(sprintf(
    "https://lpdaacsvc.cr.usgs.gov/services/gedifinder?%s=%s&%s=%s&%s=%f,%f,%f,%f&output=json",
    "version",version,
    "product", product,
    "bbox", ul_lat, ul_lon, lr_lat,lr_lon))
  content = suppressWarnings(readLines(response))
  close(response)
  results = jsonlite::parse_json(content)
  if(results$message != "") {
    stop(results$message)
  }

  out<-simplify2array(results$data)

  if (!is.null(daterange)){
    dates<-as.Date(gsub(".*(\\d{4})\\.(\\d{2})\\.(\\d{2}).*", "\\1-\\2-\\3",out))
    out_sub<-out[dates >= as.Date(daterange[1]) & dates <= as.Date(daterange[2])]

    if(length(out_sub)<1) {
      stop(paste("There are not GEDI data avalaible between",daterange[1],"and",daterange[2]))
    } else {
      return(out_sub)
    }
  } else {
    return(out)
  }
}

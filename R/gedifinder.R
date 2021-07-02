concept_ids <- list(
  GEDI01_B.001 = "C1656765475-LPDAAC_ECS",
  GEDI02_A.001 = "C1656766463-LPDAAC_ECS",
  GEDI02_B.001 = "C1656767133-LPDAAC_ECS",
  GEDI01_B.002 = "C1908344278-LPDAAC_ECS",
  GEDI02_A.002 = "C1908348134-LPDAAC_ECS",
  GEDI02_B.002 = "C1908350066-LPDAAC_ECS"
)

#' GEDI finder
#'
#' @description This function finds the exact granule(s) that contain GEDI data
#' for a given region of interest and date range
#'
#' @param product GEDI data level; Options: "GEDI01_B", "GEDI02_A" or "GEDI02_B"
#' @param ul_lat Numeric. Upper left (ul) corner coordinates, in lat
#' (decimal degrees) for the bounding box of the area of interest.
#' @param ul_lon Numeric. Upper left (ul) corner coordinates, in lon
#' (decimal degrees) for the bounding box of the area of interest.
#' @param lr_lat Numeric. Lower right (ul) corner coordinates, in lat
#' (decimal degrees) for the bounding box of the area of interest.
#' @param lr_lon Numeric. Lower right (ul) corner coordinates, in lon
#' (decimal degrees) for the bounding box of the area of interest.
#' @param version Character. The version of the GEDI product files to be
#' returned. Default "002".
#' @param daterange Vector. Date range. Specify your start and end dates
#' using ISO 8601 [YYYY]-[MM]-[DD]T[hh]:[mm]:[ss]Z. Ex.:
#' c("2019-07-01T00:00:00Z","2020-05-22T23:59:59Z"). If NULL (default),
#' the date range filter will be not applied.
#'
#' @return Return a vector object pointing out the path saving the downloaded
#' GEDI data within the boundary box coordinates provided
#'
#' @seealso bbox: Defined by the upper left and lower right corner coordinates,
#' in lat,lon ordering, for the bounding box of the area of interest
#' (e.g. \[ul_lat,ul_lon,lr_lat,lr_lon\]).
#'
#' This function relies on the existing CMR tool:
#' https://cmr.earthdata.nasa.gov/search/site/docs/search/api.html
#'
#' @examples
#' \donttest{
#' # gedifinder is a web service provided by NASA
#' # usually the request takes more than 5 seconds
#'
#' # Specifying bounding box coordinates
#' ul_lat <- 42.0
#' ul_lon <- -100
#' lr_lat <- 40.0
#' lr_lon <- -96.0
#'
#' # Specifying the date range
#' daterange <- c("2019-07-01", "2020-05-22")
#'
#' # Extracting the path to GEDI data for the specified boundary box coordinates
#' gedi02b_list <- gedifinder(
#'   product = "GEDI02_B",
#'   ul_lat,
#'   ul_lon,
#'   lr_lat,
#'   lr_lon,
#'   version = "001",
#'   daterange = daterange
#' )
#' }
#' @import jsonlite curl
#' @export
gedifinder <- function(product,
                       ul_lat,
                       ul_lon,
                       lr_lat,
                       lr_lon,
                       version = "002",
                       daterange = NULL) {
  page <- 1
  bbox <- paste(ul_lon, lr_lat, lr_lon, ul_lat, sep = ",")

  # Granules search url pattern
  url_format <- paste0(
    "https://cmr.earthdata.nasa.gov/search/granules.json?",
    "pretty=true&provider=LPDAAC_ECS&page_size=2000&concept_id=%s",
    "&bounding_box=%s"
  )
  request_url <- sprintf(
    url_format,
    concept_ids[paste0(product, ".", version)],
    bbox
  )

  # Add temporal search if not null
  if (!is.null(daterange)) {
    url_format <- paste0(request_url, "&temporal=%s,%s")
    request_url <- sprintf(url_format, daterange[1], daterange[2])
  }

  granules_href <- c()
  # Append fetched granules to granules_href
  # recursively, for each page (max 2000 per page)
  repeat {
    response <- curl::curl_fetch_memory(paste0(
      request_url,
      "&pageNum=",
      page
    ))
    content <- rawToChar(response$content)
    result <- jsonlite::parse_json(content)
    if (response$status_code != 200) {
      stop(paste("\n", result$errors, collapse = "\n"))
    }
    granules <- result$feed$entry

    if (length(granules) == 0) break

    granules_href <- c(
      granules_href,
      sapply(granules, function(x) x$links[[1]]$href)
    )
    page <- page + 1
  }

  return(granules_href)
}
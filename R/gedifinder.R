.gedi_products <- data.frame(
  product = c("GEDI01_B", "GEDI02_A", "GEDI02_B", "GEDI03", "GEDI04_A", "GEDI04_B"),
  current_version = c("003", "003", "003", "003", "003", "002.1"),
  short_name = c(
    "GEDI01_B", "GEDI02_A", "GEDI02_B",
    "GEDI_L3_LandSurface_Metrics_V3_2525",
    "GEDI_L4A_AGB_Density_V3_2508",
    "GEDI_L4B_Gridded_Biomass_V2_1_2299"
  ), stringsAsFactors = FALSE
)

.gedi_short_name <- function(product, version) {
  row <- .gedi_products[.gedi_products$product == product, , drop = FALSE]
  if (!nrow(row)) stop("Unsupported GEDI product: ", product, call. = FALSE)
  if (is.null(version)) version <- row$current_version
  normalized <- gsub("^V", "", as.character(version), ignore.case = TRUE)
  normalized <- sub("^([0-9])$", "00\\1", normalized)
  normalized <- sub("^([0-9]{2})$", "0\\1", normalized)
  if (product == "GEDI03") {
    if (normalized %in% c("002", "2")) return(list(short_name = "GEDI_L3_LandSurface_Metrics_V2_1952", version = "2"))
    return(list(short_name = row$short_name, version = "3"))
  }
  if (product == "GEDI04_A") {
    if (normalized %in% c("002.1", "2.1", "002")) return(list(short_name = "GEDI_L4A_AGB_Density_V2_1_2056", version = "2.1"))
    return(list(short_name = row$short_name, version = "3"))
  }
  if (product == "GEDI04_B") return(list(short_name = row$short_name, version = "2.1"))
  list(short_name = row$short_name, version = normalized)
}

.cmr_get <- function(url) {
  response <- curl::curl_fetch_memory(url)
  content <- jsonlite::fromJSON(rawToChar(response$content), simplifyVector = FALSE)
  if (response$status_code >= 300L) stop(paste(unlist(content$errors), collapse = "\n"), call. = FALSE)
  content
}

.gedi_link_type <- function(url) {
  if (grepl("^s3://", url)) return("s3")
  if (grepl("opendap", url, ignore.case = TRUE)) return("opendap")
  "https"
}

.select_gedi_link <- function(links, access) {
  href <- vapply(links, function(x) if (is.null(x$href)) NA_character_ else x$href, character(1))
  href <- href[!is.na(href)]
  href <- href[!grepl("(\\.xml|\\.sha256|\\.json)$", href, ignore.case = TRUE)]
  if (access == "all") return(href)
  types <- vapply(href, .gedi_link_type, character(1))
  candidates <- href[types == access]
  if (access == "https") {
    data_like <- grepl("\\.(h5|hdf5|tif|tiff|csv|zip)(\\?|$)", candidates, ignore.case = TRUE)
    if (any(data_like)) candidates <- candidates[data_like]
  }
  if (length(candidates)) candidates[[1L]] else NA_character_
}

#' Find GEDI granules through NASA CMR
#'
#' @param product One of `GEDI01_B`, `GEDI02_A`, `GEDI02_B`, `GEDI03`,
#'   `GEDI04_A`, or `GEDI04_B`.
#' @param ul_lat,ul_lon,lr_lat,lr_lon Bounding coordinates in decimal degrees.
#' @param version Product version. `NULL` selects the current supported version.
#' @param daterange Optional two-element date or date-time vector.
#' @param access Requested link type: `"https"`, `"s3"`, `"opendap"`, or `"all"`.
#' @param cloud_hosted Logical; restrict collection discovery to Earthdata Cloud.
#' @param return Return a URL vector or a metadata table.
#' @param page_size CMR page size.
#' @return A character vector or [data.table::data.table].
#' @seealso \url{https://cmr.earthdata.nasa.gov/search/site/docs/search/api.html}
#' @export
gedifinder <- function(product, ul_lat, ul_lon, lr_lat, lr_lon,
                       version = NULL, daterange = NULL,
                       access = c("https", "s3", "opendap", "all"),
                       cloud_hosted = TRUE, return = c("url", "table"),
                       page_size = 2000L) {
  access <- match.arg(access)
  return <- match.arg(return)
  product <- toupper(product)
  info <- .gedi_short_name(product, version)
  bbox <- paste(ul_lon, lr_lat, lr_lon, ul_lat, sep = ",")
  collection_url <- paste0(
    "https://cmr.earthdata.nasa.gov/search/collections.json?short_name=",
    curl::curl_escape(info$short_name), "&version=", curl::curl_escape(info$version),
    "&cloud_hosted=", tolower(as.character(isTRUE(cloud_hosted))))
  collections <- .cmr_get(collection_url)$feed$entry
  if (!length(collections)) stop("No matching GEDI collection was found in CMR.", call. = FALSE)
  collection_ids <- vapply(collections, `[[`, character(1), "id")
  rows <- list()
  k <- 0L
  for (collection_id in collection_ids) {
    page <- 1L
    repeat {
      url <- paste0(
        "https://cmr.earthdata.nasa.gov/search/granules.json?pretty=false&page_size=",
        as.integer(page_size), "&page_num=", page,
        "&collection_concept_id=", curl::curl_escape(collection_id),
        "&bounding_box=", curl::curl_escape(bbox))
      if (!is.null(daterange)) {
        if (length(daterange) != 2L) stop("'daterange' must have two values.")
        url <- paste0(url, "&temporal=", curl::curl_escape(paste(daterange, collapse = ",")))
      }
      entries <- .cmr_get(url)$feed$entry
      if (!length(entries)) break
      for (entry in entries) {
        selected <- .select_gedi_link(entry$links, access)
        if (!length(selected)) next
        for (href in selected) {
          if (is.na(href)) next
          k <- k + 1L
          rows[[k]] <- data.table::data.table(
            product = product, version = info$version, collection_id = collection_id,
            granule_id = if (is.null(entry$producer_granule_id)) basename(href) else entry$producer_granule_id,
            time_start = if (is.null(entry$time_start)) NA_character_ else entry$time_start,
            time_end = if (is.null(entry$time_end)) NA_character_ else entry$time_end,
            access = .gedi_link_type(href), url = href)
        }
      }
      if (length(entries) < page_size) break
      page <- page + 1L
    }
  }
  ans <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  if (return == "url") ans$url else ans
}

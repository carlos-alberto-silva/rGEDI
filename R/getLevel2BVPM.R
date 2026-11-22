# var.map
var.map = list()
var.map[["algorithmrun_flag"]]       = "algorithmrun_flag"
var.map[["ancillary"]]               = "ancillary"
var.map[["beam"]]                    = "beam"
var.map[["channel"]]                 = "channel"
var.map[["cover"]]                   = "cover"
var.map[["delta_time"]]              = "delta_time"
var.map[["fhd_normal"]]              = "fhd_normal"
var.map[["degrade_flag"]]            = "geolocation/degrade_flag"
var.map[["digital_elevation_model"]] = "geolocation/digital_elevation_model"
var.map[["elev_highestreturn"]]      = "geolocation/elev_highestreturn"
var.map[["elev_lowestmode"]]         = "geolocation/elev_lowestmode"
var.map[["elevation_bin0"]]          = "geolocation/elevation_bin0"
var.map[["elevation_bin0_error"]]    = "geolocation/elevation_bin0_error"
var.map[["elevation_lastbin"]]       = "geolocation/elevation_lastbin"
var.map[["elevation_lastbin_error"]] = "geolocation/elevation_lastbin_error"
var.map[["height_bin0"]]             = "geolocation/height_bin0"
var.map[["height_lastbin"]]          = "geolocation/height_lastbin"
var.map[["lat_highestreturn"]]       = "geolocation/lat_highestreturn"
var.map[["lat_lowestmode"]]          = "geolocation/lat_lowestmode"
var.map[["latitude_bin0"]]           = "geolocation/latitude_bin0"
var.map[["latitude_bin0_error"]]     = "geolocation/latitude_bin0_error"
var.map[["latitude_lastbin"]]        = "geolocation/latitude_lastbin"
var.map[["latitude_lastbin_error"]]  = "geolocation/latitude_lastbin_error"
var.map[["local_beam_azimuth"]]      = "geolocation/local_beam_azimuth"
var.map[["local_beam_elevation"]]    = "geolocation/local_beam_elevation"
var.map[["lon_highestreturn"]]       = "geolocation/lon_highestreturn"
var.map[["lon_lowestmode"]]          = "geolocation/lon_lowestmode"
var.map[["longitude_bin0"]]          = "geolocation/longitude_bin0"
var.map[["longitude_bin0_error"]]    = "geolocation/longitude_bin0_error"
var.map[["longitude_lastbin"]]       = "geolocation/longitude_lastbin"
var.map[["longitude_lastbin_error"]] = "geolocation/longitude_lastbin_error"
var.map[["solar_azimuth"]]           = "geolocation/solar_azimuth"
var.map[["solar_elevation"]]         = "geolocation/solar_elevation"
var.map[["l2a_quality_flag"]]        = "l2a_quality_flag"
var.map[["l2b_quality_flag"]]        = "l2b_quality_flag"
var.map[["land_cover_data"]]         = "land_cover_data"
var.map[["landsat_treecover"]]       = "land_cover_data/landsat_treecover"
var.map[["modis_nonvegetated"]]      = "land_cover_data/modis_nonvegetated"
var.map[["modis_nonvegetated_sd"]]   = "land_cover_data/modis_nonvegetated_sd"
var.map[["modis_treecover"]]         = "land_cover_data/modis_treecover"
var.map[["modis_treecover_sd"]]      = "land_cover_data/modis_treecover_sd"
var.map[["master_frac"]]             = "master_frac"
var.map[["master_int"]]              = "master_int"
var.map[["num_detectedmodes"]]       = "num_detectedmodes"
var.map[["omega"]]                   = "omega"
var.map[["pai"]]                     = "pai"
var.map[["pgap_theta"]]              = "pgap_theta"
var.map[["pgap_theta_error"]]        = "pgap_theta_error"
var.map[["rg"]]                      = "rg"
var.map[["rh100"]]                   = "rh100"
var.map[["rhog"]]                    = "rhog"
var.map[["rhog_error"]]              = "rhog_error"
var.map[["rhov"]]                    = "rhov"
var.map[["rhov_error"]]              = "rhov_error"
var.map[["rossg"]]                   = "rossg"
var.map[["rv"]]                      = "rv"
var.map[["selected_l2a_algorithm"]]  = "selected_l2a_algorithm"
var.map[["selected_rg_algorithm"]]   = "selected_rg_algorithm"
var.map[["sensitivity"]]             = "sensitivity"
var.map[["shot_number"]]             = "shot_number"
var.map[["stale_return_flag"]]       = "stale_return_flag"
var.map[["surface_flag"]]            = "surface_flag"


#'Get GEDI Canopy Cover and Vertical Profile Metrics (GEDI Level2B)
#'
#'@description This function extracts information from GEDI Level2B data:
#'Total Plant Area Index,	Foliage Height Diversity, Foliage Clumping Index,
#'Total Gap Probability (theta), and Total canopy cover.
#'
#'
#'@param level2b A GEDI Level2B object (output of [readLevel2B()]
#'function). An S4 object of class "gedi.level2b".
#'@param cols A character vector containing the list of columns to be extracted. See the default columns in the description.
#'@param beams Optional beam names. `NULL` reads every beam.
#'@param quality Logical; for gediDB, retain Release 3 quality observations.
#'@param collect Logical; collect remote footprints into R or retain a lazy
#'GEE/gediDB object for rasterization.
#'@param max_features Maximum GEE footprints allowed when collecting.
#'@param chunk_size Number of GEE features transferred per request.
#'@param max_granules Maximum vector granules merged by GEE.
#'
#'@return A [data.table::data.table] containing the Vegetation Profile
#'Biophysical Variables, or a lazy remote query when `collect = FALSE`.
#'
#'@seealso \url{https://www.earthdata.nasa.gov/data/catalog/lpcloud-gedi02-b-002}
#'
#'@details These are the biophysical variables and additional information extracted by default:
#'\itemize{
#'\item \emph{beam} Beam identifier
#'\item \emph{shot_number} Shot number
#'\item \emph{algorithmrun_flag} The L2B algorithm is run if this flag is set to 1 indicating data have sufficient waveform fidelity for L2B to run
#'\item \emph{l2b_quality_flag} L2B quality flag
#'\item \emph{delta_time} Transmit time of the shot since Jan 1 00:00 2018
#'\item \emph{sensitivity} Maxmimum canopy cover that can be penetrated
#'\item \emph{solar_elevation} Solar elevation
#'\item \emph{latitude_lastbin} Latitude of last bin of the pgap_theta_z, interpolated from L1B waveform coordinate
#'\item \emph{latitude_bin0} Latitude of first bin of the pgap_theta_z, interpolated from L1B waveform coordinate
#'\item \emph{elev_highestreturn} Elevation of highest detected return relative to reference ellipsoid
#'\item \emph{elev_lowestmode} Elevation of center of lowest mode relative to reference ellipsoid
#'\item \emph{rh100} RH100 slice
#'\item \emph{pai} Total Plant Area Index
#'\item \emph{fhd_normal} Foliage Height Diversity
#'\item \emph{omega} Foliage Clumping Index
#'\item \emph{pgap_theta} Total Gap Probability (theta)
#'\item \emph{cover} Total canopy cover
#'}
#'
#' Every other columns in the GEDI2B product are also available, you can specify each column by using the `cols` parameter.
#'
#'@examples
#'# Specifying the path to GEDI level2B data (zip file)
#'outdir = tempdir()
#'level2B_fp_zip <- system.file("extdata",
#'                   "GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level2A data
#'level2Bpath <- unzip(level2B_fp_zip,exdir = outdir)
#'
#'# Reading GEDI level2B data (h5 file)
#'level2b<-readLevel2B(level2Bpath=level2Bpath)
#'
#'# Extracting GEDI Vegetation Profile Biophysical Variables
#'level2BVPM<-getLevel2BVPM(level2b, cols=c("beam", "shot_number"))
#'head(level2BVPM)
#'
#'close(level2b)
#' @param parallel Logical; process a multi-granule collection concurrently.
#' @param workers Number of parallel workers, or `NULL` for automatic selection.
#'@export
getLevel2BVPM<-function(level2b, cols=c(
  "beam",
  "shot_number",
  "algorithmrun_flag",
  "l2b_quality_flag",
  "delta_time",
  "sensitivity",
  "solar_elevation",
  "latitude_lastbin",
  "latitude_bin0",
  "longitude_bin0",
  "longitude_lastbin",
  "elev_highestreturn",
  "elev_lowestmode",
  "rh100",
  "pai",
  "fhd_normal",
  "omega",
  "pgap_theta",
  "cover"
), beams = NULL, quality = FALSE, collect = TRUE,
   max_features = 100000L, chunk_size = 5000L, max_granules = 200L,
   parallel = FALSE, workers = NULL){
  if (inherits(level2b, "gedi_granule_collection")) {
    return(.extract_gedi_collection(
      level2b,
      list(cols = cols, beams = beams, quality = quality, collect = collect,
           max_features = max_features, chunk_size = chunk_size,
           max_granules = max_granules),
      parallel, workers
    ))
  }
  if (inherits(level2b, "gedi_gee_object")) {
    if (!is.null(beams)) warning("'beams' is ignored for a Google Earth Engine query; filter the returned table instead.")
    variables <- unique(as.character(cols))
    variables[variables == "l2b_quality_flag"] <- "quality_flag"
    if (!isTRUE(collect)) {
      return(.gee_select(level2b, variables, quality = quality, max_granules = max_granules))
    }
    return(.gee_collect(
      level2b, variables, quality = quality, max_features = max_features,
      chunk_size = chunk_size, max_granules = max_granules
    ))
  }
  if (inherits(level2b, "gedidb_level2b")) {
    if (!is.null(beams)) warning("'beams' is ignored for a gediDB query; filter the returned table instead.")
    requested <- unique(as.character(cols))
    aliases <- c(
      beam = "beam_name", algorithmrun_flag = "l2_algrunflag",
      l2b_quality_flag = "l2b_quality_flag_rel3", rh100 = "rh:100"
    )
    mapped <- unname(ifelse(requested %in% names(aliases), aliases[requested], requested))
    unavailable <- requested %in% c(
      "delta_time", "latitude_lastbin", "latitude_bin0",
      "longitude_bin0", "longitude_lastbin", "elev_highestreturn"
    )
    mapped <- mapped[!unavailable]
    filters <- if (isTRUE(quality)) list(l2b_quality_flag_rel3 = "== 1") else list()
    if (!isTRUE(collect)) return(.new_gedidb_query(level2b, mapped, filters))
    ans <- .gedidb_query(level2b, mapped, filters)
    restore <- c(
      beam_name = "beam", l2_algrunflag = "algorithmrun_flag",
      l2b_quality_flag_rel3 = "l2b_quality_flag"
    )
    present <- intersect(names(restore), names(ans))
    if (length(present)) data.table::setnames(ans, present, unname(restore[present]))
    return(ans[])
  }
  level2b<-level2b@h5
  groups_id<-grep("BEAM\\d{4}$",gsub("/","",
                                     .gedi_list_groups(level2b, recursive = FALSE)), value = T)
  if (!is.null(beams)) groups_id <- intersect(groups_id, beams)
  if (!length(groups_id)) return(data.table::data.table())
  m.dt<-data.table::data.table()
  pb <- utils::txtProgressBar(min = 0, max = length(groups_id), style = 3)
  i.s=0

  #i = groups_id[1]
  for ( i in groups_id){
    i.s<-i.s+1
    utils::setTxtProgressBar(pb, i.s)
    level2b_i<-level2b[[i]]
    m<-data.table::data.table()
    #col = cols[1]
    for (col in cols) {
      h5.address = var.map[[col]]
      if (is.null(h5.address)) {
        if (col %in% names(level2b_i)) {
          h5.address = col
        } else {
          if (i.s == 1) warning(sprintf("The column '%s' is not available in the GEDI2B product for BEAM %s!", col, i))
          next
        }
      }
      m[,eval(col):=level2b_i[[h5.address]][]]
    }
    m.dt<-rbind(m.dt,m)
  }
  close(pb)
  return(m.dt)
}

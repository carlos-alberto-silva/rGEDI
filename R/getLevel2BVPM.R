# l2b.var.map
l2b.var.map = list()
l2b.var.map[["algorithmrun_flag"]] = "algorithmrun_flag"
l2b.var.map[["beam"]] = "beam"
l2b.var.map[["channel"]] = "channel"
l2b.var.map[["cover"]] = "cover"
l2b.var.map[["delta_time"]] = "delta_time"
l2b.var.map[["fhd_normal"]] = "fhd_normal"
l2b.var.map[["degrade_flag"]] = "geolocation/degrade_flag"
l2b.var.map[["digital_elevation_model"]] = "geolocation/digital_elevation_model"
l2b.var.map[["elev_highestreturn"]] = "geolocation/elev_highestreturn"
l2b.var.map[["elev_lowestmode"]] = "geolocation/elev_lowestmode"
l2b.var.map[["elevation_bin0"]] = "geolocation/elevation_bin0"
l2b.var.map[["elevation_bin0_error"]] = "geolocation/elevation_bin0_error"
l2b.var.map[["elevation_lastbin"]] = "geolocation/elevation_lastbin"
l2b.var.map[["elevation_lastbin_error"]] = "geolocation/elevation_lastbin_error"
l2b.var.map[["height_bin0"]] = "geolocation/height_bin0"
l2b.var.map[["height_lastbin"]] = "geolocation/height_lastbin"
l2b.var.map[["lat_highestreturn"]] = "geolocation/lat_highestreturn"
l2b.var.map[["lat_lowestmode"]] = "geolocation/lat_lowestmode"
l2b.var.map[["latitude_bin0"]] = "geolocation/latitude_bin0"
l2b.var.map[["latitude_bin0_error"]] = "geolocation/latitude_bin0_error"
l2b.var.map[["latitude_lastbin"]] = "geolocation/latitude_lastbin"
l2b.var.map[["latitude_lastbin_error"]] = "geolocation/latitude_lastbin_error"
l2b.var.map[["local_beam_azimuth"]] = "geolocation/local_beam_azimuth"
l2b.var.map[["local_beam_elevation"]] = "geolocation/local_beam_elevation"
l2b.var.map[["lon_highestreturn"]] = "geolocation/lon_highestreturn"
l2b.var.map[["lon_lowestmode"]] = "geolocation/lon_lowestmode"
l2b.var.map[["longitude_bin0"]] = "geolocation/longitude_bin0"
l2b.var.map[["longitude_bin0_error"]] = "geolocation/longitude_bin0_error"
l2b.var.map[["longitude_lastbin"]] = "geolocation/longitude_lastbin"
l2b.var.map[["longitude_lastbin_error"]] = "geolocation/longitude_lastbin_error"
l2b.var.map[["solar_azimuth"]] = "geolocation/solar_azimuth"
l2b.var.map[["solar_elevation"]] = "geolocation/solar_elevation"
l2b.var.map[["l2a_quality_flag"]] = "l2a_quality_flag"
l2b.var.map[["l2b_quality_flag"]] = "l2b_quality_flag"
l2b.var.map[["land_cover_data"]] = "land_cover_data"
l2b.var.map[["landsat_treecover"]] = "land_cover_data/landsat_treecover"
l2b.var.map[["modis_nonvegetated"]] = "land_cover_data/modis_nonvegetated"
l2b.var.map[["modis_nonvegetated_sd"]] = "land_cover_data/modis_nonvegetated_sd"
l2b.var.map[["modis_treecover"]] = "land_cover_data/modis_treecover"
l2b.var.map[["modis_treecover_sd"]] = "land_cover_data/modis_treecover_sd"
l2b.var.map[["master_frac"]] = "master_frac"
l2b.var.map[["master_int"]] = "master_int"
l2b.var.map[["num_detectedmodes"]] = "num_detectedmodes"
l2b.var.map[["omega"]] = "omega"
l2b.var.map[["pai"]] = "pai"
l2b.var.map[["pgap_theta"]] = "pgap_theta"
l2b.var.map[["pgap_theta_error"]] = "pgap_theta_error"
l2b.var.map[["rg"]] = "rg"
l2b.var.map[["rh100"]] = "rh100"
l2b.var.map[["rhog"]] = "rhog"
l2b.var.map[["rhog_error"]] = "rhog_error"
l2b.var.map[["rhov"]] = "rhov"
l2b.var.map[["rhov_error"]] = "rhov_error"
l2b.var.map[["rossg"]] = "rossg"
l2b.var.map[["rv"]] = "rv"
l2b.var.map[["selected_l2a_algorithm"]] = "selected_l2a_algorithm"
l2b.var.map[["selected_rg_algorithm"]] = "selected_rg_algorithm"
l2b.var.map[["sensitivity"]] = "sensitivity"
l2b.var.map[["shot_number"]] = "shot_number"
l2b.var.map[["stale_return_flag"]] = "stale_return_flag"
l2b.var.map[["surface_flag"]] = "surface_flag"


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
#'
#'@return Returns an S4 object of class [data.table::data.table]
#'containing the Vegetation Profile Biophysical Variables.
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi02_bv001/
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
#'\item \emph{longitude_lastbin} Longitude of last bin of the pgap_theta_z, interpolated from L1B waveform coordinate
#'\item \emph{longitude_bin0} Longitude of first bin of the pgap_theta_z, interpolated from L1B waveform coordinate
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
#'A special column `pavd` is available, which is the sum of the pavd_z.
#'
#'The following columns are also available:
#'\itemize{
#'\item \emph{channel} Channel number (0-7)
#'\item \emph{master_frac} Master time, fractional part. master_int+master_frac is equivalent to /BEAMXXXX/geolocation/delta_time and /BEAMXXXX/geophys_corr/delta_time.
#'\item \emph{master_int} Master time, integer part. Seconds since master_time_epoch. master_int+master_frac is equivalent to /BEAMXXXX/geolocation/delta_time and /BEAMXXXX/geophys_corr/delta_time.
#'\item \emph{num_detectedmodes} Number of detected modes in rxwaveform
#'\item \emph{pgap_theta_error} Error of the estimated Pgap(theta) for the selected L2A algorithm
#'\item \emph{l2a_quality_flag} Flag simplifying selection of most useful data
#'\item \emph{rg} Integral of the ground component in the RX waveform for the selected L2A processing version
#'\item \emph{rhog} Volumetric scattering coefficient of the ground (reflectance x phase function)
#'\item \emph{rhog_error} Error term in Rho (ground)
#'\item \emph{rhov} Volumetric scattering coefficient of the canopy (reflectance x phase function)
#'\item \emph{rhov_error} Error term in Rho (canopy)
#'\item \emph{rossg} Mean projection of unit leaf area on a plane perpendicular to the direction of the laser beam at view zenith angle theta
#'\item \emph{rv} Integral of the vegetation component in the RX waveform for the selected L2A processing version
#'\item \emph{selected_l2a_algorithm} 	ID of algorithm selected as identifying the lowest non-noise mode
#'\item \emph{selected_rg_algorithm} 	0 = L2B algorithm not run; 1 = algorithm 1 (successful); 2 = algorithm 1 (partial success - valid tx_eg parameters were unavailable); 3 = algorithm 2
#'\item \emph{stale_return_flag} Indicates that a "stale" cue point from the coarse search algorithm is being used.
#'\item \emph{surface_flag} Indicates elev_lowestmode is within 300m of DEM or MSS
#'\item \emph{degrade_flag} Greater than zero if the shot occurs during a degrade period, zero otherwise.
#'\item \emph{digital_elevation_model} Digital elevation model height above the WGS84 ellipsoid. Interpolated at latitude_bin0 and longitude_bin0 from the TandemX 90m product.
#'\item \emph{elevation_bin0} Elevation of first bin of the pgap_theta_z, interpolated from L1B waveform coordinate
#'\item \emph{elevation_bin0_error} Error in elevation of first bin of the pgap_theta_z, interpolated from L1B waveform coordinate
#'\item \emph{elevation_lastbin} Elevation of last bin of the pgap_theta_z, interpolated from L1B waveform coordinate
#'\item \emph{elevation_lastbin_error} Error in elevation of last bin of the pgap_theta_z, interpolated from L1B waveform coordinate
#'\item \emph{height_bin0} Height of the first bin of the pgap_theta_z, relative to the ground.
#'\item \emph{height_lastbin} Height of the last bin of the pgap_theta_z, relative to the ground.
#'\item \emph{lat_highestreturn} Latitude of highest detected return
#'\item \emph{lat_lowestmode} Latitude of center of lowest mode
#'\item \emph{latitude_bin0_error} Error in latitude of first bin of the pgap_theta_z, interpolated from L1B waveform coordinate
#'\item \emph{latitude_lastbin_error} Error in latitude of last bin of the pgap_theta_z, interpolated from L1B waveform coordinate
#'\item \emph{local_beam_azimuth} Azimuth of the unit pointing vector for the laser in the local ENU frame. The angle is measured from North and positive towards East.
#'\item \emph{local_beam_elevation} Elevation of the unit pointing vector for the laser in the local ENU frame. The angle is measured from East-North plane and positive towards Up.
#'\item \emph{lon_highestreturn} Longitude of highest detected return
#'\item \emph{lon_lowestmode} Longitude of center of lowest mode
#'\item \emph{longitude_bin0_error} Error in longitude of first bin of the pgap_theta_z, interpolated from L1B waveform coordinate
#'\item \emph{longitude_lastbin_error} Error in longitude of last bin of the pgap_theta_z, interpolated from L1B waveform coordinate
#'\item \emph{solar_azimuth} The azimuth of the sun position vector from the laser bounce point position in the local ENU frame. The angle is measured from North and is positive towards East.
#'\item \emph{landsat_treecover} Tree cover in the year 2010, defined as canopy closure for all vegetation taller than 5m in height (Hansen et al.). Encoded as a percentage per output grid cell.
#'\item \emph{modis_nonvegetated} Percent non-vegetated from MODIS data. Interpolated at latitude_bin0 and longitude_bin0. doi:10.5067/MODIS/MOD44B.006
#'\item \emph{modis_nonvegetated_sd} Percent non-vegetated standard deviation from MODIS data. Interpolated at latitude_bin0 and longitude_bin0. doi:10.5067/MODIS/MOD44B.006
#'\item \emph{modis_treecover} Percent tree cover from MODIS data. Interpolated at latitude_bin0 and longitude_bin0. doi:10.5067/MODIS/MOD44B.006
#'\item \emph{modis_treecover_sd} Percent tree cover standard deviation from MODIS data. Interpolated at latitude_bin0 and longitude_bin0. doi:10.5067/MODIS/MOD44B.006
#'}
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
#'@export
getLevel2BVPM <- function(level2b, cols = c(
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
)) {
  level2b <- level2b@h5
  groups_id <- grep("BEAM\\d{4}$", gsub("/", "",
                                     hdf5r::list.groups(level2b, recursive = F)), value = T)
  m.dt <- data.table::data.table()
  pb <- utils::txtProgressBar(min = 0, max = length(groups_id), style = 3)
  i_s = 0

  for (i in groups_id) {
    i_s = i_s + 1
    utils::setTxtProgressBar(pb, i.s)
    level2b_i = level2b[[i]]
    m = data.table::data.table()
    for (col in cols) {
      h5_address = l2b.var.map[[col]]
      if (col == "pavd")  {
        h5_address = "pavd_z"
        m[, eval(col) := colSums(level2b_i[[h5_address]][, ])]
        next
      }
      if (is.null(h5_address)) {
        if (level2b_i$exists(col)) {
        h5_address = col
        } else {
          if (i.s == 1) warning(
            sprintf(
              "The column '%s' is not available in the GEDI2B product!",
              col
            )
          )
          m[, eval(col) := NA]
          next
        }
      }
      base_addr = gsub("^(.*)/.*", "\\1", h5_address)
      if (level2b_i$exists(base_addr) && level2b_i$exists(h5_address))
        m[, eval(col) := level2b_i[[h5_address]][]]
    }
    m_dt = data.table::rbindlist(list(m.dt, m), fill = TRUE)
  }
  close(pb)
  return(m_dt)
}

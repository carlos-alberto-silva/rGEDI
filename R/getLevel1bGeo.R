#'Get Pulse FullWaveform Geolocation (GEDI Level1B)
#'
#'@description This function extracts Pulse FullWaveform Geolocation from GEDI Level1B data
#'
#'@usage getLevel1BGeo(level1b, select = c("latitude_bin0", "latitude_lastbin", "longitude_bin0", "longitude_lastbin", "shot_number") )
#'
#'@param level1b A GEDI Level1B object (output of \code{\link[rGEDI:readLevel1B]{readLevel1B}} function). A S4 object of class "gedi.level1b".
#'@param select A character vector specifying the fields to extract from GEDI Level1B data. The default is set to extract \emph{latitude_bin0}, \emph{latitude_lastbin}, \emph{longitude_bin0}, \emph{longitude_lastbin}, and \emph{shot_number}. See details for more options.
#'
#'@return A S4 object of class \code{\link[data.table:data.table]{data.table-class}}.
#'
#'@details Besides of the fields specified as defaul, GEDI mission level 1B description (https://lpdaac.usgs.gov/products/gedi01_bv001/) details all the information that can be extracted:
#'\itemize{
#'\item \emph{all_samples_sum	Sum} - of all values within the 10 km range window.
#'\item \emph{beam} -	Beam number	Number.
#'\item \emph{channel} -	Channel number.
#'\item \emph{master_frac} -	Master time, fractional part.
#'\item \emph{master_int} -	Master time, integer part.
#'\item \emph{noise_mean_corrected} -	Noise mean.
#'\item \emph{noise_stddev_corrected} -	Corrected noise standard deviation.
#'\item \emph{nsemean_even} -	Noise mean of the beam's detector channel from even sub-converte.
#'\item \emph{nsemean_odd} -	Noise mean of the beam's odd sub-converter.
#'\item \emph{rx_energy} -	Integrated energy in receive (RX) waveform after subtracting the noise mean.
#'\item \emph{rx_offset} -	Time interval from first stored sample to first downloaded RX sample.
#'\item \emph{rx_open} -	Time interval from time 0 to first stored RX sample.
#'\item \emph{rx_sample_count} -	The number of sample intervals (elements) in each RX waveform.
#'\item \emph{rx_sample_start_index} -	The index in the rxwaveform dataset of the first element of each RX waveform starting at 1.
#'\item \emph{rxwaveform} -	The corrected receive RX waveforms using rx_sample_count and rx_sample_start_index to identify the location of each waveform.
#'\item \emph{selection_stretchers_x} -	Commanded number of samples added to algorithm section on the left.
#'\item \emph{selection_stretchers_y} -	Commanded number of samples added to algorithm section on the right.
#'\item \emph{shot_number} -	Unique shot identifier.
#'\item \emph{stale_return_flag} -	Indicates that a "stale" cue point from the coarse search algorithm is being used.
#'\item \emph{th_left_used} -	Count values for left threshold used in fine search where two consecutive points at or above this value indicates pulse detection.
#'\item \emph{tx_egamplitude} -	Amplitude of the extended Gaussian fit to the transmit (TX) waveform.
#'\item \emph{tx_egamplitude_error} -	Error on tx_egamplitude.
#'\item \emph{tx_egbias} -	Bias of the extended Gaussian fit to the TX waveform.
#'\item \emph{tx_egbias_error} -	Error on tx_egbias.
#'\item \emph{tx_egflag} -	Extended Gaussian fit status flag.
#'\item \emph{tx_eggamma} -	Gamma value of the extended Gaussian fit to the TX waveform.
#'\item \emph{tx_eggamma_error} -	Error on tx_eggamma.
#'\item \emph{tx_egsigma} -	Sigma of the extended Gaussian fit to the TX waveform.
#'\item \emph{tx_egsigma_error} -	Error on tx_egsigma.
#'\item \emph{tx_gloc} -	Location (mean) of the Gaussian fit to the TX waveform.
#'\item \emph{tx_gloc_error} -	Error on tx_gloc.
#'\item \emph{tx_pulseflag} -	Set to 1 if a pulse is detected in the TX waveform.
#'\item \emph{tx_sample_count} -	The number of sample intervals (elements) in each transmit waveform.
#'\item \emph{tx_sample_start_index} -	The index in the rxwaveform dataset of the first element of each RX waveform starting at 1.
#'\item \emph{txwaveform} -	The corrected TX waveforms using tx_sample_count and tx_sample_start_index to identify the location of each waveform.
#'\item \emph{master_time_epoch} -	Number of GPS seconds between the GPS epoch (1980-01-06T00:00:00Z) and the GEDI epoch (2018-01-01T00:00:00Z).
#'\item \emph{mean_samples} -	Averaging window (shots) for non-repeated odd-even measurements from housekeeping.
#'\item \emph{smoothing_width} -	The 1-sigma width of the gaussian kernel used to smooth the corrected waveform.
#'\item \emph{altitude_instrument} -	Height of the instrument diffractive optical element (DOE) above the WGS84 ellipsoid.
#'\item \emph{altitude_instrument_error} -	Error on altitude_instrument.
#'\item \emph{bounce_time_offset_bin0} -	The difference between the TX time and the time at the start of the RX window.
#'\item \emph{bounce_time_offset_bin0_error} -	Error on bounce_time_offset_bin0.
#'\item \emph{bounce_time_offset_lastbin} -	The difference between the TX time and the time at the end of the RX window.
#'\item \emph{bounce_time_offset_lastbin_error} -	Error on bounce_time_offset_lastbin.
#'\item \emph{degrade} -	Greater than zero if the shot occurs during a degrade period, zero otherwise.
#'\item \emph{delta_time} -	Transmit time of the shot, measured in seconds since 2018-01-01.
#'\item \emph{digital_elevation_model} -	Digital elevation model height above the WGS84 ellipsoid.
#'\item \emph{elevation_bin0} -	Height of the start of the RX window, relative to the WGS-84 ellipsoid.
#'\item \emph{elevation_bin0_error} -	Error on elevation_bin0.
#'\item \emph{elevation_lastbin} -	Height of the end of the RX window, relative to the WGS-84 ellipsoid.
#'\item \emph{elevation_lastbin_error} -	Error on elevation_lastbin.
#'\item \emph{latitude_bin0} -	Latitude of the start of the RX window.
#'\item \emph{latitude_bin0_error} -	Error on latitude_bin0.
#'\item \emph{latitude_lastbin} -	Latitude of the end of the RX window.
#'\item \emph{latitude_lastbin_error} -	Error on latitude_lastbin.
#'\item \emph{latitude_instrument} -	Latitude of the instrument diffractive optical element (DOE) at laser transmit time.
#'\item \emph{latitude_instrument_error} -	Error on latitude_instrument.
#'\item \emph{local_beam_azimuth} -	Azimuth of the unit pointing vector for the laser in the local East, North, Up (ENU) frame.
#'\item \emph{local_beam_azimuth_error} -	Error on local_beam_azimuth.
#'\item \emph{local_beam_elevation} -	Elevation of the unit pointing vector for the laser in the local ENU frame.
#'\item \emph{local_beam_elevation_error} -	Error on local_beam_elevation.
#'\item \emph{longitude_bin0} -	Longitude of the start of the RX window.
#'\item \emph{longitude_bin0_error} -	Error on longitude_bin0.
#'\item \emph{longitude_lastbin} -	Longitude of the end of the RX window.
#'\item \emph{longitude_lastbin_error} -	Error on longitude_lastbin.
#'\item \emph{longitude_instrument} -	Longitude of the instrument diffractive optical element (DOE) at laser transmit time.
#'\item \emph{longitude_instrument_error} -	Error on longitude_instrument.
#'\item \emph{mean_sea_surface} -	Mean sea surface height above the WGS84 ellipsoid, includes the geoid	.
#'\item \emph{neutat_delay_derivative_bin0} -	Change in neutral atmospheric delay per height change for the start of the RX window.
#'\item \emph{neutat_delay_derivative_lastbin} -	Change in neutral atmospheric delay per height change for the end of the RX window.
#'\item \emph{neutat_delay_total_bin0} -	Total neutral atmosphere delay correction (wet+dry) from the TX pulse tothe start of the RX window.
#'\item \emph{neutat_delay_total_lastbin} -	Total neutral atmosphere delay correction (wet+dry) from the TX pulse tothe end of the RX window.
#'\item \emph{range_bias_correction} -	The range bias applied to the range measurement.
#'\item \emph{shot_number} -	Unique shot identifier	Number.
#'\item \emph{solar_azimuth} -	The azimuth of the sun position vector.
#'\item \emph{solar_elevation} -	The elevation of the sun position vector.
#'\item \emph{surface_type} -	Flags describing which surface types.
#'\item \emph{delta_time} -	Transmit time of the shot.
#'\item \emph{dynamic_atmosphere_correction} -	Dynamic Atmospheric Correction (DAC) includes inverted barometer (IB) effect.
#'\item \emph{geoid} -	Geoid height above WGS-84 reference ellipsoid.
#'\item \emph{tide_earth} -	Solid Earth tides.
#'\item \emph{tide_load} -	Load Tide - Local displacement due to Ocean Loading.
#'\item \emph{tide_ocean} -	Ocean Tides including diurnal and semi-diurnal, and longerperiod tides.
#'\item \emph{tide_ocean_pole} -	Oceanic surface rotational deformation due to polar motion.
#'\item \emph{tide_pole} -	Solid Earth Pole Tide -- Rotational deformation due to polar motion.
#'}
#'
#'@examples
#'# specify the path to GEDI Level 1B data
#'level1bpath <- system.file("extdata", "GEDIexample_level01B.h5", package="rGEDI")
#'
#'# Reading GEDI level1B data
#'level1b <- readLevel1B(level1bpath)
#'
#'# Get GEDI level1B geolocations
#'level1bGeo<-getLevel1BGeo(level1b,select=c("latitude_bin0", "longitude_bin0","shot_number"))
#'head(level1bGeo)
#'
#'
#'@export
getLevel1BGeo<-function(level1b,select=c("latitude_bin0", "latitude_lastbin", "longitude_bin0", "longitude_lastbin", "shot_number")) {
  level1b<-level1b@h5

  datasets<-hdf5r::list.datasets(level1b, recursive = T)
  datasets_names<-basename(datasets)

  selected<-datasets_names %in% select

  for ( i in select){
    if  ( i =="shot_number"){
      assign(i,bit64::as.integer64(NaN))
    } else {
      assign(i,numeric())
    }
  }

  dtse2<-datasets[selected][!grepl("geolocation/shot_number",datasets[selected])]

  # Set progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(dtse2), style = 3)
  i.s=0


  # i = "BEAM0000/shot_number"
  for ( i in dtse2){
    #print(i)
    i.s<-i.s+1
    utils::setTxtProgressBar(pb, i.s)
    name_i<-basename(i)
    if ( name_i =="shot_number"){

      assign(name_i, bit64::c.integer64(get(name_i),level1b[[i]][]))

    } else {

      assign(name_i, c(get(name_i), level1b[[i]][]))

    }

  }

  if ( "shot_number" %in% select){
    level1b.dt<-data.table::data.table(get("shot_number"))[-1,]
    colnames(level1b.dt)<-"shot_number"
    select2<-select[!select[]=="shot_number"]

  } else{
    level1b.dt<-data.table::data.table(get(select[1]))
    colnames(level1b.dt)<-select[1]
    select2<-select[-1]
  }

  for ( i in select2){
    level1b.dt[,i]<-get(i)
  }

  close(pb)
  return(level1b.dt)
}


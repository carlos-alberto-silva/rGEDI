#'Clip GEDI Level1B data by Coordinates
#'
#'@description This function clips GEDI Level1B data (geolocated waveforms) within a given bounding coordinates
#'
#'@usage clipLevel1B(level1b, xmin, xmax, ymin, ymax, output)
#'
#'@param level1b A [`gedi.level1b-class`] object (output of [readLevel1B()] function).
#'An S4 object of class [`gedi.level1b-class`].
#'@param xmin Numeric. West longitude (x) coordinate of the bounding rectangle, in decimal degrees.
#'@param xmax Numeric. East longitude (x) coordinate of the bounding rectangle, in decimal degrees.
#'@param ymin Numeric. South latitude (y) coordinate of the bounding rectangle, in decimal degrees.
#'@param ymax Numeric. North latitude (y) coordinate of the bounding rectangle, in decimal degrees.
#'@param output Optional character path where to save the new hdf5file. The default stores a temporary file only.
#'
#'@return Returns a list of S4 objects of class [`gedi.level1b-class`] containing clipped GEDI Level1B data.
#'
#'@seealso \url{https://lpdaac.usgs.gov/products/gedi01_bv002/}
#'
#'@examples
#'\donttest{
#'# Specifying the path to GEDI level1B data (zip file)
#'outdir = tempdir()
#'
#'level1B_fp_zip <- system.file("extdata",
#'                   "GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level1B data
#'level1Bpath <- unzip(level1B_fp_zip,exdir = outdir)
#'
#'# Reading GEDI level1B data (h5 file)
#'level1b<-readLevel1B(level1Bpath=level1Bpath)
#'
#'# Bounding rectangle coordinates
#'xmin=-44.13
#'xmax=-44.12
#'ymin=-13.74
#'ymax=-13.73
#'
#'# Specifying output file and path
#'output<-file.path(outdir,"GEDI01_B_2019108080338_O01964_T05337_02_003_01_clip")
#'
#'# Clipping GEDI Level1B data by extent boundary box
#'level1b_clip <- clipLevel1B(level1b,xmin, xmax, ymin, ymax,output)
#'
#' close(level1b)
#' close(level1b_clip)
#'}
#'@import hdf5r fs
#'@export
clipLevel1B = function(level1b, xmin, xmax, ymin, ymax, output=""){
  output = checkOutput(output)
  checkClipExtentInputs(level1b, "gedi.level1b", xmin, xmax, ymin, ymax)

  # Get all spatial data as a list of dataframes with spatial information
  spData = getSpatialData1B(level1b)

  masks = clipSpDataByExtentLevelB(spData, xmin, xmax, ymin, ymax)

  newFile = clipByMask1B(level1b,
                         masks,
                         output)
  output = newFile@h5$filename
  close(newFile)
  result = readLevel1B(output)

  return (result)
}

#'Clip GEDI Level1B data by geometry
#'
#'@description This function clips GEDI Level1B (geolocated waveforms) data within a given bounding geometry
#'
#'
#'@param level1b A [`gedi.level1b-class`] object (output of [readLevel1B()] function).
#'An S4 object of class "gedi.level1b".
#'@param polygon_spdf Polygon. An object of class [`sp::SpatialPolygonsDataFrame`][sp::SpatialPolygonsDataFrame-class] from `sp` package,
#'which can be loaded as an ESRI shapefile using [raster::shapefile()].
#'@param output Optional character path where to save the new [`hdf5r::H5File`][hdf5r::H5File-class]. The default stores a temporary file only.
#'@param split_by Polygon id. If defined, GEDI data will be clipped by each polygon using the attribute specified by `split_by` from the attribute table.
#'
#'@return Returns a list of S4 object of class [`gedi.level1b-class`] containing clipped GEDI Level1B data.
#'
#'@examples
#'\donttest{
#'outdir = tempdir()
#'
#'# Specifying the path to GEDI level1B data (zip file)
#'level1B_fp_zip <- system.file("extdata",
#'                   "GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level1B data
#'level1Bpath <- unzip(level1B_fp_zip,exdir = outdir)
#'
#'# Reading GEDI level1B data (h5 file)
#'level1b<-readLevel1B(level1Bpath=level1Bpath)
#'
#'# Specifying the path to shapefile
#'polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package="rGEDI")
#'
#'# Reading shapefile as SpatialPolygonsDataFrame object
#'library(raster)
#'polygon_spdf<-shapefile(polygon_filepath)
#'
#'# Spepecifing output file and path
#'output<-file.path(outdir,"GEDI01_B_2019108080338_O01964_T05337_02_003_01_clip")
#'
#'# clipping GEDI Level1B data by extent boundary box
#'level1b_clip <- clipLevel1BGeometry(level1b, polygon_spdf = polygon_spdf,
#'                                    output=output,
#'                                    split_by="id")
#'
#'close(level1b)
#'lapply(level1b_clip, close)
#'}
#'@import hdf5r
#'@export
clipLevel1BGeometry = function(level1b, polygon_spdf, output="", split_by=NULL) {
  output = checkOutput(output)
  checkClipGeoInputs(level1b, "gedi.level1b", polygon_spdf, split_by)

  spData = getSpatialData1B(level1b)

  xmin = polygon_spdf@bbox[1,1]
  xmax = polygon_spdf@bbox[1,2]
  ymin = polygon_spdf@bbox[2,1]
  ymax = polygon_spdf@bbox[2,2]

  masks = clipSpDataByExtentLevelB(spData, xmin, xmax, ymin, ymax)
  polygon_masks = getPolygonMaskLevelB(spData, masks, polygon_spdf, split_by)
  results = clipByMasks(level1b, polygon_masks, output, split_by, clipByMask1B)

  return (results)
}


# Helper function to return spatial data within a dataframe
getSpatialData1B = function(level1b) {
  level1b.h5<-level1b@h5
  groups_id<-grep("BEAM\\d{4}$",gsub("/","",
                                     hdf5r::list.groups(level1b.h5, recursive = F)), value = T)

  beams_spdf = list()

  for ( i in gsub("/","",groups_id)){
    beams_spdf[[i]] = data.frame(
      latitude_bin0=level1b.h5[[paste0(i,"/geolocation/latitude_bin0")]][],
      latitude_lastbin=level1b.h5[[paste0(i,"/geolocation/latitude_lastbin")]][],
      longitude_bin0=level1b.h5[[paste0(i,"/geolocation/longitude_bin0")]][],
      longitude_lastbin=level1b.h5[[paste0(i,"/geolocation/longitude_lastbin")]][]
    )
  }

  return (beams_spdf)
}

clipByMask1B = function(level1b, masks, output = "") {
  newFile =  hdf5r::H5File$new(output, mode="w")
  if(length(masks) == 0) {
    message("No intersection found!")
    newFile$close_all()
    newFile = hdf5r::H5File$new(output, mode="r")
    result = new("gedi.level1b", h5 = newFile)
    return(result)
  }

  for (attr in hdf5r::list.attributes(level1b@h5)) {
    hdf5r::h5attr(newFile, attr) = hdf5r::h5attr(level1b@h5, attr)
  }

  all_groups = hdf5r::list.groups(level1b@h5)
  beams_with_value = sapply(masks, length)>0
  beams_with_value = names(which(beams_with_value))
  beams_with_value = c(beams_with_value, "METADATA")
  which_groups = gsub("([^/]*).*","\\1",all_groups) %in% beams_with_value
  groups_with_value = all_groups[which_groups]

  # Setup progress bar
  all_datasets = hdf5r::list.datasets(level1b@h5)
  which_datasets = gsub("([^/]*).*","\\1",all_datasets) %in% beams_with_value
  datasets_with_value = all_datasets[which_datasets]

  total = length(datasets_with_value)
  pb = utils::txtProgressBar(min = 0, max = total, style = 3)
  progress = 0

  for (group in groups_with_value) {
    beam_id = strsplit(group, "/")[[1]][1]
    mask = masks[[beam_id]]
    mask_size = length(mask)

    hdf5r::createGroup(newFile, group)
    createAttributesWithinGroup(level1b@h5, newFile, group)

    if(beam_id != "METADATA") {
      beam_shot_n = level1b@h5[[beam_id]][["shot_number"]]$dims
      total_waveforms = list()
      total_waveforms[["rx"]] = sum(level1b@h5[[beam_id]][["rx_sample_count"]][])
      total_waveforms[["tx"]] = sum(level1b@h5[[beam_id]][["tx_sample_count"]][])
    }


    for (dt in hdf5r::list.datasets(level1b@h5[[group]], recursive = FALSE, full.names = T)) {
      if (grepl("[rt]x_sample_start_index$", dt)) {
        progress = progress + 1
        utils::setTxtProgressBar(pb, progress)
        next
      }
      h5_dt = level1b@h5[[dt]]
      dt_dim = h5_dt$dims
      dtype = h5_dt$get_type()
      if (is.na(all(h5_dt$chunk_dims))) {
        chunkdims = NULL
      } else {
        chunkdims = h5_dt$chunk_dims
      }


      if (length(dt_dim) == 1) {
        if (dt_dim == 1) {
          hdf5r::createDataSet(newFile,dt,h5_dt[], dtype=dtype, chunk_dim=chunkdims)
        } else if (dt_dim == beam_shot_n) {
          hdf5r::createDataSet(newFile,dt,h5_dt[mask], dtype=dtype, chunk_dim=chunkdims)
        } else if (dt_dim %in% total_waveforms | dt_dim %% beam_shot_n == 0) {
          prefix = ifelse(substr(basename(dt),1,2)=="rx", "rx", "tx")
          sampleCount = sprintf("%s_sample_count", prefix)
          sampleStartIndex = sprintf("%s_sample_start_index", prefix)
          countPath=file.path(beam_id, sampleCount, fsep = "/")
          startIdxPath=file.path(beam_id, sampleStartIndex, fsep = "/")

          counts = level1b@h5[[countPath]][mask]
          startIdx = level1b@h5[[startIdxPath]][mask]
          v.seq = Vectorize(seq.default,vectorize.args = c("from", "length.out"), SIMPLIFY=T)
          mask_list_waveform=v.seq(startIdx, length.out=counts)
          mask_waveform=Reduce(c, mask_list_waveform)
          newStartIdx = c(1,cumsum(counts)+1)[-mask_size-1]
          newFile$create_dataset(
            name       = startIdxPath,
            robj       = newStartIdx,
            dims       = length(newStartIdx),
            chunk_dims = level1b@h5[[startIdxPath]]$chunk_dims,
            dtype      = level1b@h5[[startIdxPath]]$get_type())

          total_size = length(mask_waveform)
          chunk_part = 1
          dt_res=hdf5r::createDataSet(newFile, dt, dtype=dtype, chunk_dim=chunkdims, dims=total_size)
          while (chunk_part < total_size) {
            end = chunk_part+chunkdims-1
            if (end > total_size) {
              end = total_size
            }
            get_part = mask_waveform[(chunk_part):(end)]
            dt_res[chunk_part:end] =  h5_dt[get_part]
            chunk_part = end+1
          }
        }
      } else if (length(dt_dim) == 2 && dt_dim[1] == beam_shot_n) {
        if (length(mask) == 1) {
          chunkdims = chunkdims[[2]]
        }
        hdf5r::createDataSet(newFile,dt,level1b@h5[[dt]][mask,], dtype=dtype, chunk_dim=chunkdims)
      } else {
        stop(paste0("Don't know how to handle dataset: ", dt, "\nContact the maintainer of the package!"))
      }

      #Update progress
      progress = progress + 1
      utils::setTxtProgressBar(pb, progress)
    }
  }

  newFile$close_all()
  newFile =  hdf5r::H5File$new(output, mode="r")
  result = new("gedi.level1b", h5 = newFile)
  close(pb)
  #spatial = level1B2dt(level1b)
  return (result)
}

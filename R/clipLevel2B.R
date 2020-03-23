#'Clip GEDI Level2B data by Coordinates
#'
#'@description This function extracts GEDI Level1B data a within given bounding coordinates
#'
#'
#'@param level2b A GEDI Level2B object (output of \code{\link[rGEDI:readLevel2B]{readLevel2B}} function).
#'An S4 object of class "gedi.level2b".
#'@param xmin Numeric. West longitude (x) coordinate of the bounding rectangle, in decimal degrees.
#'@param xmax Numeric. East longitude (x) coordinate of the bounding rectangle, in decimal degrees.
#'@param ymin Numeric. South latitude (y) coordinate of the bounding rectangle, in decimal degrees.
#'@param ymax Numeric. North latitude (y) coordinate of the bounding rectangle, in decimal degrees.
#'@param output Optional character path where to save the new hdf5 file. The default stores a temporary file only.
#'
#'@return Returns a list of S4 object of class "gedi.level2b" containing clipped GEDI Level2B data.
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi01_bv001/
#'
#'@examples
#'\donttest{
#'outdir = tempdir()
#'
#'# Specifying the path to GEDI level2B data (zip file)
#'level2B_fp_zip <- system.file("extdata",
#'                   "GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level2A data
#'level2Bpath <- unzip(level2B_fp_zip,exdir = dirname(level2B_fp_zip))
#'
#'# Reading GEDI level2B data (h5 file)
#'level2b<-readLevel2B(level2Bpath=level2Bpath)
#'
#'# Bounding rectangle coordinates
#'xmin=-44.13
#'xmax=-44.12
#'ymin=-13.74
#'ymax=-13.73
#'
#'# Specifying output file and path
#'output<-file.path(outdir,"GEDI02_B_2019108080338_O01964_T05337_02_001_01_clip")
#'
#'# Clipping GEDI data by extent boundary box
#'level2b_clip <- clipLevel2B(level2b, xmin, xmax, ymin, ymax)
#'
#'close(level2b)
#'close(level2b_clip)
#'}
#'@export
clipLevel2B = function(level2b, xmin, xmax, ymin, ymax, output=""){
  # Get all spatial data as a list of dataframes with spatial information
  spData = getSpatialData2B(level2b)
  checkClipExtentInputs(level2b, "gedi.level2b", xmin, xmax, ymin, ymax)

  masks = clipSpDataByExtentLevelB(spData, xmin, xmax, ymin, ymax)

  if (output == "") {
    output = tempfile(fileext = ".h5")
  }
  output = fs::path_ext_set(output, "h5")

  newFile = clipByMask2B(level2b,
                         masks,
                         output)
  output = newFile@h5$filename
  close(newFile)
  result = readLevel2B(output)

  return (result)
}

#'Clip GEDI Level2B data by geometry
#'
#'@description This function extracts GEDI Level1B data within a given geometry
#'
#'@param level2b A GEDI Level2B object (output of \code{\link[rGEDI:readLevel2B]{readLevel2B}} function).
#'An S4 object of class "gedi.level2b".
#'@param polygon_spdf Polygon. An object of class \code{\link[sp]{SpatialPolygonsDataFrame-class}},
#'which can be loaded as an ESRI shapefile using \code{\link[raster:shapefile]{raster::shapefile()}} function in the \emph{raster} package.
#'@param output optional character path where to save the new h5file. Default "" (temporary file).
#'@param split_by Polygon id. If defined, GEDI data will be clipped by each polygon using the attribute specified by \code{split_by} from the attribute table.
#'
#'@return Returns a list of S4 objects of class "gedi.level2b" containing clipped GEDI Level2B data.
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi01_bv001/
#'
#'@examples
#'\donttest{
#'outdir = tempdir()
#'
#'# Specifying the path to GEDI level2B data (zip file)
#'level2B_fp_zip <- system.file("extdata",
#'                   "GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level2A data
#'level2Bpath <- unzip(level2B_fp_zip,exdir = dirname(level2B_fp_zip))
#'
#'# Reading GEDI level2B data (h5 file)
#'level2b<-readLevel2B(level2Bpath=level2Bpath)
#'
#'# Specifying the path to shapefile
#'polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package="rGEDI")
#'
#'# Reading shapefile as SpatialPolygonsDataFrame object
#'library(raster)
#'polygon_spdf<-shapefile(polygon_filepath)
#'
#'# Specifying output file and path
#'output<-file.path(outdir,"GEDI02_B_2019108080338_O01964_T05337_02_001_01_clip")
#'
#'# Clipping GEDI data by extent boundary box
#'level2b_clip <- clipLevel2BGeometry(level2b, polygon_spdf = polygon_spdf,
#'                                    output=output,
#'                                    split_by="id")
#'
#'close(level2b)
#'lapply(level2b_clip, close)
#'}
#'@export
clipLevel2BGeometry = function(level2b, polygon_spdf, output="", split_by=NULL) {
  output = checkOutput(output)
  checkClipGeoInputs(level2b, "gedi.level2b", polygon_spdf, split_by)

  spData = getSpatialData2B(level2b)

  xmin = polygon_spdf@bbox[1,1]
  xmax = polygon_spdf@bbox[1,2]
  ymin = polygon_spdf@bbox[2,1]
  ymax = polygon_spdf@bbox[2,2]

  masks = clipSpDataByExtentLevelB(spData, xmin, xmax, ymin, ymax)
  polygon_masks = getPolygonMaskLevelB(spData, masks, polygon_spdf, split_by)
  results = clipByMasks(level2b, polygon_masks, output, split_by, clipByMask2B)

  return (results)
}


# Helper function to return spatial data within a dataframe
getSpatialData2B = function(level2b) {
  level2b.h5<-level2b@h5
  groups_id<-grep("BEAM\\d{4}$",gsub("/","",
                                     hdf5r::list.groups(level2b.h5, recursive = F)), value = T)

  beams_spdf = list()

  for ( i in gsub("/","",groups_id)){
    beams_spdf[[i]] = data.frame(
      latitude_bin0=level2b.h5[[paste0(i,"/geolocation/latitude_bin0")]][],
      latitude_lastbin=level2b.h5[[paste0(i,"/geolocation/latitude_lastbin")]][],
      longitude_bin0=level2b.h5[[paste0(i,"/geolocation/longitude_bin0")]][],
      longitude_lastbin=level2b.h5[[paste0(i,"/geolocation/longitude_lastbin")]][]
    )
  }
  return (beams_spdf)
}


clipByMask2B = function(level2b, masks, output = "") {
  newFile =  hdf5r::H5File$new(output, mode="w")
  if(length(masks) == 0) {
    message("No intersection found!")
    newFile$close_all()
    newFile = hdf5r::H5File$new(output, mode="r")
    result = new("gedi.level2b", h5 = newFile)
    return(result)
  }

  for (attr in hdf5r::list.attributes(level2b@h5)) {
    hdf5r::h5attr(newFile, attr) = hdf5r::h5attr(level2b@h5, attr)
  }

  all_groups = hdf5r::list.groups(level2b@h5)
  beams_with_value = sapply(masks, length)>0
  beams_with_value = names(which(beams_with_value))
  beams_with_value = c(beams_with_value, "METADATA")
  which_groups = gsub("([^/]*).*","\\1",all_groups) %in% beams_with_value
  groups_with_value = all_groups[which_groups]

  # Setup progress bar
  all_datasets = hdf5r::list.datasets(level2b@h5)
  which_datasets = gsub("([^/]*).*","\\1",all_datasets) %in% beams_with_value
  datasets_with_value = all_datasets[which_datasets]

  total = length(datasets_with_value)
  pb = utils::txtProgressBar(min = 0, max = total, style = 3)
  progress = 0

  for (group in groups_with_value) {
    beam_id = strsplit(group, "/")[[1]][1]
    mask = masks[[beam_id]]
    mask_size = length(mask)

    hdf5r::createGroup(newFile,group)
    createAttributesWithinGroup(level2b@h5, newFile, group)

    for (dt in hdf5r::list.datasets(level2b@h5[[group]], recursive = FALSE, full.names = T)) {
      beam_shot_n = level2b@h5[[beam_id]][["shot_number"]]$dims
      h5_dt = level2b@h5[[dt]]
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
        } else if ((dt_dim %% beam_shot_n) == 0) {
          n_waveforms = h5_dt$dims / beam_shot_n
          v.seq = Vectorize(seq.default,vectorize.args = c("from"), SIMPLIFY=T)
          mask_init = mask*n_waveforms - (n_waveforms - 1)
          mask_waveform = matrix(v.seq(mask_init, len=n_waveforms), nrow=1)[1,]
          total_size = n_waveforms*mask_size
          chunk_part = 1
          dt_res=hdf5r::createDataSet(newFile, dt, dtype=dtype, chunk_dim=chunkdims, dims=total_size)
          while (chunk_part < total_size) {
            end = chunk_part+chunkdims-1
            if (end > total_size) {
              end = total_size
            }
            get_part = mask_waveform[(chunk_part):(end)]
            dt_res[get_part] =  h5_dt[get_part]
            chunk_part = end+1
          }
        }
      } else if (length(dt_dim) == 2 && dt_dim[1] == beam_shot_n) {
        if (length(mask) == 1) {
          chunkdims = chunkdims[[2]]
        }
        hdf5r::createDataSet(newFile,dt,h5_dt[mask,], dtype=dtype, chunk_dim=chunkdims)
      } else if (length(dt_dim) == 2 && dt_dim[2] == beam_shot_n){
        if (length(mask) == 1) {
          chunkdims = chunkdims[[1]]
        }
        hdf5r::createDataSet(newFile,dt,h5_dt[,mask], dtype=dtype, chunk_dim=chunkdims)
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
  result = new("gedi.level2b", h5 = newFile)
  close(pb)
  #spatial = level2b2dt(level2b)
  return (result)
}

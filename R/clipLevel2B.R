#'Clip GEDI Level2B data by Coordinates
#'
#'@description This function extracts GEDI Level1B data within given bounding coordinates
#'
#'
#'@param level2b A GEDI Level2B object (output of \code{\link[rGEDI:readLevel2B]{readLevel2B}} function). A S4 object of class "gedi.level2b".
#'@param xmin Numeric. West longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param xmax Numeric. East longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param ymin Numeric. South latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'@param ymax Numeric. North latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'@param output Optional character path where to save the new hdf5 file. The default stores a temporary file only.
#'@return An S4 object of class "gedi.level2b".
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi01_bv001/
#'
#'@examples
#'# specify the path and data file and read it
#'level2bpath <- system.file("extdata", "GEDIexample_level02B.h5", package="rGEDI")
#'
#'# reading GEDI level2B data
#'level2b <- readLevel2B(level2bpath)
#'
#'# Bounding rectangle coordinates
#'xmin = -44.15036
#'xmax = -44.10066
#'ymin = -13.75831
#'ymax = -13.71244
#'
#'# clip level2BVPM by extent boundary box
#'level2b_clip <- level2BVPM(level2BVPM,xmin, xmax, ymin, ymax)
#'
#'library(leaflet)
#'leaflet() %>%
#'  addCircleMarkers(level2b_clip@dt$longitude_bin0,
#'                   level2b_clip@dt$latitude_bin0,
#'                   radius = 1,
#'                   opacity = 1,
#'                   color = "red")  %>%
#'  addScaleBar(options = list(imperial = FALSE)) %>%
#'  addPolygons(data=polygon_spdf,weight=1,col = 'white',
#'              opacity = 1, fillOpacity = 0) %>%
#'  addProviderTiles(providers$Esri.WorldImagery)
#'@export
clipLevel2B = function(level2b, xmin, xmax, ymin, ymax, output=""){
  # Get all spatial data as a list of dataframes with spatial information
  spData = getSpatialData2B(level2b)

  masks = clipSpDataByExtentLevelB(spData, xmin, xmax, ymin, ymax)

  if (output == "") {
    output = tempfile(fileext = ".h5")
  }
  output = fs::path_ext_set(output, "h5")

  newFile = clipByMask2B(level2b,
                         masks,
                         output)
  output = newFile@h5$filename
  newFile@h5$close_all()
  result = readLevel2B(output)

  return (result)
}

#'Clip GEDI Level2B data by geometry
#'
#'@description This function extracts GEDI Level1B data within given geometry
#'
#'@param level2b h5file; S4 object of class H5File
#'@param polygon_spdf Polygon. An object of class \code{\link[sp]{SpatialPolygonsDataFrame-class}},
#'which can be loaded as an ESRI shapefile using \code{\link[rgdal:readOGR]{readOGR}} function in the \emph{rgdal} package.
#'@param output optional character path where to save the new h5file. Default "" (temporary file).
#'@return An S4 object of class "gedi.level2b".
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi01_bv001/
#'
#'@examples
#'
#'# specify the path and data file and read it
#'level2bpath <- system.file("extdata", "GEDIexample_level02B.h5", package="rGEDI")
#'
#'# reading GEDI level2B data
#'level2b <- readLevel2B(level2bpath)
#'
#'# specify the path to shapefile
#'polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package="rGEDI")
#'
#'# Reading shapefile as SpatialPolygonsDataFrame object
#'library(rgdal)
#'polygon_spdf<-readOGR(polygons_filepath)
#'
#'# clip level2BVPM by geometry
#'level2b_clip_geometry <- clipLevel2BGeometry(level2BVPM,polygon_spdf=polygon_spdf)
#'
#'@export
clipLevel2BGeometry = function(level2b, polygon_spdf, output="", split_by=NULL) {
  output = checkOutput(output)
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

  for (attr in hdf5r::list.attributes(level2b@h5)) {
    hdf5r::h5attr(newFile, attr) = hdf5r::h5attr(level2b@h5, attr)
  }

  all_groups = hdf5r::list.groups(level2b@h5)

  # Setup progress bar
  all_datasets = hdf5r::list.datasets(level2b@h5)
  total = length(all_datasets)
  pb = utils::txtProgressBar(min = 0, max = total, style = 3)
  progress = 0

  for (group in all_groups) {
    beam_id = strsplit(group, "/")[[1]][1]
    mask = masks[[beam_id]]
    mask_size = length(mask)

    if (mask_size == 0)  {
      n_datasets = length(hdf5r::list.datasets(level2b@h5[[group]]))
      progress = progress + n_datasets
      utils::setTxtProgressBar(pb, progress)
      next
    }

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
        stop(paste0("Don't know how to treat dataset: ", dt, "\nContact the maintainer of the package!"))
      }

      #Update progress
      progress = progress + 1
      utils::setTxtProgressBar(pb, progress)
    }
  }

  level2b@h5 = newFile
  close(pb)
  #spatial = level2b2dt(level2b)
  return (level2b)
}

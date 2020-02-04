#'Clip GEDI Level1B data
#'
#'@description Clip GEDI Level1 data within a given bounding coordinates
#'
#'
#'@param level1b h5file; S4 object of class H5File
#'@param xleft numeric. left x coordinates of rectangles.
#'@param xright numeric. right x coordinates of rectangles.
#'@param ybottom numeric. bottom y coordinates of rectangles.
#'@param ytop numeric. top y coordinates of rectangles.
#'@param output optional character path where to save the new h5file. Default "" (temporary file).
#'
#'@return Returns An object of class H5File; subset of LVIS Level1 data
#'@author Caio Hamamura
#'@examples
#'
#'#' LVIS level 2 file path
#' #level1bpath = system.file("extdata", "lvis_level1_clip.h5", package="rLVIS")
#'
#'# Rectangle
#' #xleft = 82
#' #xright = 82.8
#' #ybottom = 0
#' #ytop = 20
#'
#'#' Reading LVIS level 2 file
#' #level1_waveform = readLevel1B(level1bpath)
#'
#' #output = tempfile(fileext=".h5")
#'
#' #clipped_waveform = clipLevel1B(level1_waveform, output, xleft, xright, ybottom, ytop)
#'
#'@export
#'
clipLevel1Bh5 = function(level1b, xleft, xright, ybottom, ytop, output=""){
  # Get all spatial data as a list of dataframes with spatial information
  spData = getSpatialData(level1b)

  masks = lapply(spData, function(x) {
    mask = x$longitude_bin0 >= xleft &
      x$longitude_bin0 <= xright &
      x$latitude_bin0 >= ybottom &
      x$latitude_bin0 <= ytop &
      x$longitude_lastbin >= xleft &
      x$longitude_lastbin <= xright &
      x$latitude_lastbin >= ybottom &
      x$latitude_lastbin <= ytop

    return ((1:length(x$longitude_bin0))[mask])
  })

  newFile = clipByMask(level1b,
                       masks,
                       output)
  output = newFile@h5$filename
  hdf5r::h5close(newFile@h5)
  result = readLevel1B(output)

  return (result)
}

#'Clip LVIS Level1 data by geometry
#'
#'@description Clip LVIS Level1 data within a given bounding coordinates
#'
#'
#'@param level1b h5file; S4 object of class H5File
#'@param polygon_spdf SpatialDataFrame. A polygon dataset for clipping the waveform
#'@param output optional character path where to save the new h5file. Default "" (temporary file).
#'
#'@return Returns An object of class H5File; subset of LVIS Level1 data
#'@author Caio Hamamura
#'@examples
#'
#'#' LVIS level 2 file path
#' #level1_filepath = system.file("extdata", "lvis_level1_clip.h5", package="rLVIS")
#'
#'#' Reading LVIS level 2 file
#' #level1_waveform = readLevel1b(level1_filepath)
#'
#'# Polgons file path
#' #polygons_filepath <- system.file("extdata", "LVIS_Mondah_clip_polygon.shp", package="rLVIS")
#'
#'# Reading LVIS level 2 file
#' #polygon_spdf<-raster::shapefile(polygons_filepath)
#'
#' #output = tempfile(fileext="h5")
#'
#' #clipped_waveform = clipLevel1Geometry(level1_waveform, output, polygon_spdf)
#'
#'@export
clipLevel1Bh5Geometry = function(level1b, polygon_spdf, output="") {
  spData = getSpatialData(level1b)

  xleft = polygon_spdf@bbox[1,1]
  xright = polygon_spdf@bbox[1,2]
  ybottom = polygon_spdf@bbox[2,1]
  ytop = polygon_spdf@bbox[2,2]

  masks = lapply(spData, function(x) {
    mask = x$longitude_bin0 >= xleft &
      x$longitude_bin0 <= xright &
      x$latitude_bin0 >= ybottom &
      x$latitude_bin0 <= ytop &
      x$longitude_lastbin >= xleft &
      x$longitude_lastbin <= xright &
      x$latitude_lastbin >= ybottom &
      x$latitude_lastbin <= ytop

    return ((1:length(x$longitude_bin0))[mask])
  })

  message("Intersecting with polygon...")
  pb = utils::txtProgressBar(min = 0, max = length(masks), style = 3)
  progress = 0
  for (beam in names(masks)) {
    mask = masks[[beam]]

    if (length(mask) == 0) next

    spDataMasked = spData[[beam]][mask,]
    points = sp::SpatialPointsDataFrame(coords=matrix(c(spDataMasked$longitude_bin0, spDataMasked$latitude_bin0), ncol=2),
                                        data=data.frame(id=mask), proj4string = polygon_spdf@proj4string)
    pts = raster::intersect(points, polygon_spdf)
    masks[[beam]] = as.integer(pts@data$id)

    progress = progress + 1
    utils::setTxtProgressBar(pb, progress)
  }
  close(pb)

  message("Writing new HDF5 file...")
  newFile = clipByMask(level1b,
                       masks,
                       output)
  output = newFile@h5$filename
  hdf5r::h5close(newFile@h5)
  result = readLevel1B(output)

  return (result)
}


# Helper function to return spatial data within a dataframe
getSpatialData = function(level1b) {
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

clipByMask = function(level1b, masks, output = "") {
  if (output == "") {
    tmp_filename <- tempfile(fileext = ".h5")
    newFile =  hdf5r::H5File$new(tmp_filename, mode="w")
  } else {
    newFile =  hdf5r::H5File$new(output, mode="w")
  }

  for (attr in hdf5r::list.attributes(level1b@h5)) {
    hdf5r::h5attr(newFile, attr) = hdf5r::h5attr(level1b@h5, attr)
  }


  all_groups = hdf5r::list.groups(level1b@h5)

  # Setup progress bar
  all_datasets = hdf5r::list.datasets(level1b@h5)
  total = length(all_datasets)
  pb = utils::txtProgressBar(min = 0, max = total, style = 3)
  progress = 0

  for (group in all_groups) {
    beam_id = strsplit(group, "/")[[1]][1]
    mask = masks[[beam_id]]
    mask_size = length(mask)

    if (mask_size == 0)  {
      n_datasets = length(hdf5r::list.datasets(level1b@h5[[group]]))
      progress = progress + n_datasets
      utils::setTxtProgressBar(pb, progress)
      next
    }

    hdf5r::createGroup(newFile,group)

    for (dt in hdf5r::list.datasets(level1b@h5[[group]], recursive = FALSE, full.names = T)) {
      beam_shot_n = level1b@h5[[beam_id]][["shot_number"]]$dims
      dt_dim = level1b@h5[[dt]]$dims

      if (length(dt_dim) == 1) {
        if (dt_dim == 1) {
          hdf5r::createDataSet(newFile,dt,level1b@h5[[dt]][])
        } else if (dt_dim == beam_shot_n) {
          hdf5r::createDataSet(newFile,dt,level1b@h5[[dt]][mask])
        } else if ((dt_dim %% beam_shot_n) == 0) {
          n_waveforms = level1b@h5[[dt]]$dims / beam_shot_n
          v.seq = Vectorize(seq.default,vectorize.args = c("from"), SIMPLIFY=T)
          mask_init = mask*n_waveforms - (n_waveforms - 1)
          mask_waveform = matrix(v.seq(mask_init, len=n_waveforms), nrow=1)[1,]
          waveform=level1b@h5[[dt]][mask_waveform]
          hdf5r::createDataSet(newFile,dt,waveform)
        }
      } else if (length(dt_dim) == 2 && dt_dim[1] == beam_shot_n) {
        hdf5r::createDataSet(newFile,dt,level1b@h5[[dt]][mask,])
      } else {
        stop(paste0("Don't know how to treat dataset: ", dt, "\nContact the maintainer of the package!"))
      }

      #Update progress
      progress = progress + 1
      utils::setTxtProgressBar(pb, progress)
    }
  }

  level1b@h5 = newFile
  #spatial = getLevel1bGeo(level1b)
  return (level1b)
}

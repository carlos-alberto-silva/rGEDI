#'Clip GEDI Level1B data
#'
#'@description Clip LVIS Level1 data within a given bounding coordinates
#'
#'
#'@param level1_waveform h5file; S4 object of class H5File
#'@param output path where to save the new h5file
#'@param xleft numeric. left x coordinates of rectangles.
#'@param xright numeric. right x coordinates of rectangles.
#'@param ybottom numeric. bottom y coordinates of rectangles.
#'@param ytop numeric. top y coordinates of rectangles.
#'
#'@return Returns An object of class H5File; subset of LVIS Level1 data
#'@author Caio Hamamura
#'@examples
#'
#'#' LVIS level 2 file path
#'level1_filepath = system.file("extdata", "lvis_level1_clip.h5", package="rLVIS")
#'
#'# Rectangle
#'xleft = 9.35986
#'xright = 9.35988
#'ybottom = 0.5786
#'ytop = 0.5790
#'
#'#' Reading LVIS level 2 file
#'level1_waveform = readLevel1b(level1_filepath)
#'
#'output = tempfile(fileext="h5")
#'
#'clipped_waveform = clipLevel1(level1_waveform, output, xleft, xright, ybottom, ytop)
#'
#'@export
#'
clipLevel1 = function(level1b, output, xleft=-120.14159, xright=-90.70369, ybottom=46.83897, ytop=54.40533){

  spData = getSpatialData(level1b)


  # xleft ybottom xright ytop
  mask =
    spData$longitude_bin0 >= xleft &
    spData$longitude_bin0 <= xright &
    spData$latitude_bin0 >= ybottom &
    spData$latitude_bin0 <= ytop &
    spData$longitude_lastbin >= xleft &
    spData$longitude_lastbin <= xright &
    spData$latitude_lastbin >= ybottom &
    spData$latitude_lastbin <= ytop

  mask = (1:length(spData$longitude_bin0))[mask]

  newFile = clipByMask(level1b,
                     #  output,
                       mask)

  return (newFile)
}

#'Clip LVIS Level1 data by geometry
#'
#'@description Clip LVIS Level1 data within a given bounding coordinates
#'
#'
#'@param level1_waveform h5file; S4 object of class H5File
#'@param output path where to save the new h5file
#'@param polygon_spdf SpatialDataFrame. A polygon dataset for clipping the waveform
#'
#'@return Returns An object of class H5File; subset of LVIS Level1 data
#'@author Caio Hamamura
#'@examples
#'
#'#' LVIS level 2 file path
#'level1_filepath = system.file("extdata", "lvis_level1_clip.h5", package="rLVIS")
#'
#'#' Reading LVIS level 2 file
#'level1_waveform = readLevel1b(level1_filepath)
#'
#'# Polgons file path
#'polygons_filepath <- system.file("extdata", "LVIS_Mondah_clip_polygon.shp", package="rLVIS")
#'
#'# Reading LVIS level 2 file
#'polygon_spdf<-raster::shapefile(polygons_filepath)
#'
#'output = tempfile(fileext="h5")
#'
#'clipped_waveform = clipLevel1Geometry(level1_waveform, output, polygon_spdf)
#'
#'@export
clipLevel1Geometry = function(level1_waveform, output, polygon_spdf) {
  spData = getSpatialData(level1b)

  points = sp::SpatialPointsDataFrame(coords=matrix(c(spData$longitude_bin0, spData$latitude_bin0), ncol=2),
                                      data=data.frame(id=1:length(spData$longitude_bin0)), proj4string = polygon_spdf@proj4string)
  pts = raster::intersect(points, polygon_spdf)
  mask = as.integer(pts@data$id)

  newFile = clipByMask(level1_waveform,
                       output,
                       mask)

  return (newFile)
}


getSpatialData = function(level1b) {

  level1b<-level1b@h5
  groups_id<-grep("BEAM\\d{4}$",gsub("/","",
                                     hdf5r::list.groups(level1b, recursive = F)), value = T)
  latitude_bin0<-NULL
  latitude_lastbin<-NULL
  longitude_bin0<-NULL
  longitude_lastbin<-NULL
  for ( i in gsub("/","",groups_id)){
    latitude_bin0=c(latitude_bin0,level1b[[paste0(i,"/geolocation/latitude_bin0")]][])
    latitude_lastbin=c(latitude_lastbin,level1b[[paste0(i,"/geolocation/latitude_lastbin")]][])
    longitude_bin0=c(longitude_bin0,level1b[[paste0(i,"/geolocation/longitude_bin0")]][])
    longitude_lastbin=c(longitude_lastbin,level1b[[paste0(i,"/geolocation/longitude_lastbin")]][])
  }

  return (list(latitude_bin0=latitude_bin0,
               longitude_bin0=longitude_bin0,
               latitude_lastbin=latitude_lastbin,
               longitude_lastbin=longitude_lastbin))
}

clipByMask = function(level1b, output, mask) {

  test_filename <- tempfile(fileext = ".h5")
  newFile =  hdf5r::H5File$new(test_filename, mode="w")

  for (attr in hdf5r::list.attributes(level1b@h5)) {
    hdf5r::h5attr(newFile, attr) = hdf5r::h5attr(level1b@h5, attr)
  }

  for (group in hdf5r::list.groups(level1b@h5)) {
    hdf5r::createGroup(newFile,group)

    for (dt in hdf5r::list.datasets(level1b@h5[[group]], recursive = FALSE, full.names = T)) {
      hdf5r::createDataSet(newFile,dt,level1b@h5[[dt]][])
    }
  }


  for (dataset in h5::list.datasets(level1_waveform, recursive = F)) {
    newFile[dataset] = level1_waveform[dataset][mask]
  }

  # Adjust ancillary extents if present
  if (h5::existsGroup(newFile, "ancillary_data")) {
    if (!h5::existsDataSet(newFile, "ancillary_data/Maximum Latitude")) {
      h5::createDataSet(newFile, "ancillary_data/Maximum Latitude", type="character", dimensions=1)
      h5::createDataSet(newFile, "ancillary_data/Maximum Longitude", type="character", dimensions=1)
      h5::createDataSet(newFile, "ancillary_data/Minimum Latitude", type="character", dimensions=1)
      h5::createDataSet(newFile, "ancillary_data/Minimum Longitude", type="character", dimensions=1)
    }
    maxY = max(newFile["LAT0"][], newFile["LAT1023"][])
    maxX = max(newFile["LON0"][], newFile["LON1023"][])
    minY = min(newFile["LAT0"][], newFile["LAT1023"][])
    minX = min(newFile["LON0"][], newFile["LON1023"][])

    writeH5field(newFile["ancillary_data/Maximum Latitude"], sprintf("%.8f",maxY))
    writeH5field(newFile["ancillary_data/Maximum Longitude"], sprintf("%.8f",maxX))
    writeH5field(newFile["ancillary_data/Minimum Latitude"], sprintf("%.8f",minY))
    writeH5field(newFile["ancillary_data/Minimum Longitude"], sprintf("%.8f",minX))
  }

  return (newFile)
}


writeH5field = function(dataset, value) {
  ds = h5::selectDataSpace(dataset, 1)
  h5::writeDataSet(dataset, value, ds)
}

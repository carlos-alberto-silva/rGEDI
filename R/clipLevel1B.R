#'Clip GEDI Level1B data by Coordinates
#'
#'@description This function clips GEDI Level1B data within given bounding coordinates
#'
#'@usage clipLevel1B(level1b, xleft, xright, ybottom, ytop, output)
#'
#'@param level1b A GEDI Level1B object (output of \code{\link[rGEDI:readLevel1B]{readLevel1B}} function). A S4 object of class "gedi.level1b".
#'@param xleft Numeric. West longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param xright Numeric. East longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param ybottom Numeric. South latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'@param ytopNumeric. North latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'@param output Optional character path where to save the new hdf5file. The default stores a temporary file only.
#'
#'@return An S4 object of class "gedi.level1b".
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi01_bv001/
#'
#'@examples
#'# specify the path to GEDI data
#'level1bpath <- system.file("extdata", "GEDIexample_level01B.h5", package="rGEDI")
#'
#'# Reading GEDI level1B file
#'level1b<-readLevel1b(level1bpath)
#'
#'# Bounding rectangle coordinates
#'xleft = -116.4683
#'xright = -116.5583
#'ybottom = 46.75208
#'ytop = 46.84229
#'
#'# clip by extent boundary box
#'level1b_clip <- clipLevel1B(level1b,xleft, xright, ybottom, ytop)
#'
#'@import hdf5r
#'@export
clipLevel1B = function(level1b, xleft, xright, ybottom, ytop, output=""){
  if (output == "") {
    output = tempfile(fileext = ".h5")
  }

  # Get all spatial data as a list of dataframes with spatial information
  spData = getSpatialData1B(level1b)

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

  newFile = clipByMask1B(level1b,
                       masks,
                       output)
  output = newFile@h5$filename
  hdf5r::h5close(newFile@h5)
  result = readLevel1B(output)

  return (result)
}

#'Clip GEDI Level1B data by geometry
#'
#'@description This function clips GEDI Level1B data within a bounding geometry
#'
#'@usage clipLevel1BGeometry(level1b, polygon_spdf, output)
#'
#'@param level1b A GEDI Level1B object (output of \code{\link[rGEDI:readLevel1B]{readLevel1B}} function). A S4 object of class "gedi.level1b".
#'@param polygon_spdf Polygon. An object of class \code{\link[sp]{SpatialPolygonsDataFrame-class}},
#'which can be loaded as an ESRI shapefile using \code{\link[rgdal:readOGR]{readOGR}} function in the \emph{rgdal} package.
#'@param output Optional character path where to save the new hdf5file. The default stores a temporary file only.
#'
#'@return An S4 object of class "gedi.level1b".
#'
#'@examples
#'# specify the path to GEDI data
#'level1bpath <- system.file("extdata", "GEDIexample_level01B.h5", package="rGEDI")
#'
#'# Reading GEDI level1B file
#'level1b<-readLevel1b(level1bpath)
#'
#'# specify the path to shapefile
#'polygon_filepath <- system.file("extdata", "clip_polygon.shp", package="rGEDI")
#'
#'# Reading shapefile as SpatialPolygonsDataFrame object
#'library(rgdal)
#'polygon_spdf<-readOGR(polygons_filepath)
#'
#'# clip by extent boundary box
#'level1b_clip <- clipLevel1BGeometry(level1b, polygon_spdf = polygon_spdf)
#'
#'@import hdf5r
#'@export
clipLevel1BGeometry = function(level1b, polygon_spdf, output="") {
  if (output == "") {
    output = tempfile(fileext = ".h5")
  }

  spData = getSpatialData1B(level1b)

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
  polygon_masks = list()

  for (beam in names(masks)) {
    mask = masks[[beam]]

    if (length(mask) == 0) next

    spDataMasked = spData[[beam]][mask,]
    points = sp::SpatialPointsDataFrame(coords=matrix(c(spDataMasked$longitude_bin0, spDataMasked$latitude_bin0), ncol=2),
                                        data=data.frame(id=mask), proj4string = polygon_spdf@proj4string)
    pts = raster::intersect(points, polygon_spdf)
    for (pol_id in levels(pts@data$d)) {
      polygon_masks[[pol_id]][[beam]] = pts[pts@data$d == pol_id,]@data[,1]
    }

    progress = progress + 1
    utils::setTxtProgressBar(pb, progress)
  }
  close(pb)

  message("Writing new HDF5 file...")
  results = list()
  for (pol_id in names(polygon_masks)) {
    output2 = gsub("\\.h5$", paste0("_", pol_id,".h5"), output)
    results[[pol_id]] = clipByMask1B(level1b,
                                     polygon_masks[[pol_id]],
                                     output2)
  }

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
  #spatial = level1B2dt(level1b)
  return (level1b)
}

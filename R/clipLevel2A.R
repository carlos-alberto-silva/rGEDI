#'Clip GEDI level2a data
#'
#'@description Clip GEDI Level1 data within a given bounding coordinates
#'
#'
#'@param level2a h5file; S4 object of class H5File
#'@param xleft numeric. left x coordinates of rectangles.
#'@param xright numeric. right x coordinates of rectangles.
#'@param ybottom numeric. bottom y coordinates of rectangles.
#'@param ytop numeric. top y coordinates of rectangles.
#'@param output optional character path where to save the new h5file. Default "" (temporary file).
#'
#'@return Returns An object of class H5File; subset of LVIS Level2 data
#'@author Caio Hamamura
#'@examples
#'
#'#' LVIS level 2A file path
#' #level2apath = system.file("extdata", "lvis_level1_clip.h5", package="rLVIS")
#'
#'# Rectangle
#' #xleft = 81
#' #xright = 83
#' #ybottom = 2
#' #ytop = 4
#'
#'#' Reading LVIS level 2 file
#' #level1_waveform = readlevel2a(level2apath)
#'
#' #output = tempfile(fileext=".h5")
#'
#' #clipped_waveform = cliplevel2a(level1_waveform, output, xleft, xright, ybottom, ytop)
#'
#'@export
#'
cliplevel2ah5 = function(level2a, xleft, xright, ybottom, ytop, output=""){
  if (output == "") {
    output = tempfile(fileext = ".h5")
  }

  # Get all spatial data as a list of dataframes with spatial information
  spData = getSpatialData2A(level2a)

  masks = lapply(spData, function(x) {
    masks2 = lapply(x, function(y) {
      mask = y$longitude_lowest >= xleft &
        y$longitude_lowest <= xright &
        y$latitude_highest >= ybottom &
        y$latitude_highest <= ytop &
        y$longitude_highest >= xleft &
        y$longitude_highest <= xright &
        y$latitude_lowest >= ybottom &
        y$latitude_lowest <= ytop

      return ((1:length(y$longitude_lowest))[mask])
    })
    return (masks2)
  })

  newFile = clipByMask2A(level2a,
                       masks,
                       output)
  output = newFile@h5$filename
  hdf5r::h5close(newFile@h5)
  result = readlevel2a(output)

  return (result)
}

#'Clip LVIS Level1 data by geometry
#'
#'@description Clip LVIS Level1 data within a given bounding coordinates
#'
#'
#'@param level2a h5file; S4 object of class H5File
#'@param polygon_spdf SpatialDataFrame. A polygon dataset for clipping the waveform
#'@param output optional character path where to save the new h5file. Default "" (temporary file).
#'
#'@return Returns An object of class H5File; subset of LVIS Level1 data
#'@author Caio Hamamura
#'@examples
#'
#'#' LVIS level 2 file path
#' #level1_filepath = system.file("extdata", "biomas.zip", package="rGEDI")
#'
#'#' Reading LVIS level 2 file
#' #level1_waveform = readlevel2a(level1_filepath)
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
cliplevel2ah5Geometry = function(level2a, polygon_spdf, output="") {

  if (output == "") {
    output = tempfile(fileext = ".h5")
  }

  spData = getSpatialData2A(level2a)

  xleft = polygon_spdf@bbox[1,1]
  xright = polygon_spdf@bbox[1,2]
  ybottom = polygon_spdf@bbox[2,1]
  ytop = polygon_spdf@bbox[2,2]

  masks = lapply(spData, function(x) {
    masks2 = lapply(x, function(y) {
      mask = y$longitude_lowest >= xleft &
        y$longitude_lowest <= xright &
        y$latitude_highest >= ybottom &
        y$latitude_highest <= ytop &
        y$longitude_highest >= xleft &
        y$longitude_highest <= xright &
        y$latitude_lowest >= ybottom &
        y$latitude_lowest <= ytop

      return ((1:length(y$longitude_lowest))[mask])
    })
    return (masks2)
  })

  message("Intersecting with polygon...")
  pb = utils::txtProgressBar(min = 0, max = length(masks), style = 3)
  progress = 0
  polygon_masks = list()

  for (beam in names(masks)) {
    masks2 = masks[[beam]]

    for (i in 1:length(masks2)) {
      mask = masks2[[i]]
      if (length(mask) == 0) next

      spDataMasked = spData[[beam]][[i]][mask,]
      points = sp::SpatialPointsDataFrame(coords=matrix(c(spDataMasked$longitude_highest, spDataMasked$latitude_highest), ncol=2),
                                          data=data.frame(id=mask), proj4string = polygon_spdf@proj4string)
      pts = raster::intersect(points, polygon_spdf)
      for (pol_id in levels(pts@data$d)) {
        mask_name = names(masks2)[i]
        polygon_masks[[pol_id]][[beam]][[mask_name]] = pts[pts@data$d == pol_id,]@data[,1]
      }

      progress = progress + 1
      utils::setTxtProgressBar(pb, progress)
    }
  }
  close(pb)

  message("Writing new HDF5 file...")

  results = list()
  for (pol_id in names(polygon_masks)) {
    output2 = gsub("\\.h5$", paste0("_", pol_id,".h5"), output)
    results[[pol_id]] = clipByMask2A(level2a,
                                 polygon_masks[[pol_id]],
                                 output2)
  }
  names(results) = NULL

  return (results)
}


# Helper function to return spatial data within a dataframe
getSpatialData2A = function(level2a) {
  level2a.h5<-level2a@h5
  groups_id<-grep("BEAM\\d{4}$",gsub("/","",
                                     hdf5r::list.groups(level2a.h5, recursive = F)), value = T)

  beams_spdf = list()

  for ( i in gsub("/","",groups_id)){
    algs_coordinates = list()


    n_algs = level2a.h5[[i]][["ancillary/l2a_alg_count"]][1]


    beams_spdf[[i]] = algs_coordinates

    for (j in 1:n_algs) {
      beams_spdf[[i]][[j]] = data.frame(
        latitude_highest=level2a.h5[[paste0(i,"/geolocation/lat_highestreturn_a", j)]][],
        latitude_lowest=level2a.h5[[paste0(i,"/geolocation/lat_lowestreturn_a", j)]][],
        longitude_highest=level2a.h5[[paste0(i,"/geolocation/lon_highestreturn_a", j)]][],
        longitude_lowest=level2a.h5[[paste0(i,"/geolocation/lon_lowestreturn_a", j)]][]
      )
    }

    beams_spdf[[i]][["main"]] = data.frame(
      latitude_highest=level2a.h5[[i]][["lat_highestreturn"]][],
      latitude_lowest=level2a.h5[[i]][["lat_lowestmode"]][],
      longitude_highest=level2a.h5[[i]][["lon_highestreturn"]][],
      longitude_lowest=level2a.h5[[i]][["lon_lowestmode"]][]
    )
  }

  return (beams_spdf)
}

clipByMask2A = function(level2a, masks, output = "") {
  newFile =  hdf5r::H5File$new(output, mode="w")

  for (attr in hdf5r::list.attributes(level2a@h5)) {
    hdf5r::h5attr(newFile, attr) = hdf5r::h5attr(level2a@h5, attr)
  }


  all_groups = hdf5r::list.groups(level2a@h5)

  # Setup progress bar
  all_datasets = hdf5r::list.datasets(level2a@h5)
  total = length(all_datasets)
  pb = utils::txtProgressBar(min = 0, max = total, style = 3)
  progress = 0

  # Check if the beam has any intersecting area
  groups_have_value = lapply(lapply(masks, function(x) sapply(x, length)), sum)>0

  for (group in all_groups) {
    beam_id = strsplit(group, "/")[[1]][1]

    if (group %in% names(groups_have_value) && ! groups_have_value[group])  {
      n_datasets = length(hdf5r::list.datasets(level2a@h5[[group]]))
      progress = progress + n_datasets
      utils::setTxtProgressBar(pb, progress)
      next
    }

    hdf5r::createGroup(newFile, group)
    createAttributesWithinGroup(level2a@h5, newFile, group)

    # Datasets to loop
    datasets = hdf5r::list.datasets(level2a@h5[[group]], recursive = FALSE, full.names = T)

    # Create list of algorithm ids for the datasets to choose right mask
    alg_ids = as.list(as.integer(sapply(regmatches(datasets, regexec("_a(\\d+)", datasets)), function(x) x[2])))
    names(alg_ids) = datasets

    for (dt in datasets) {
      # Get right mask
      alg_id = alg_ids[[dt]]

      if (is.na(alg_id)) {
        mask = masks[[beam_id]][["main"]]
      } else {
        mask = masks[[beam_id]][[alg_id]]
      }
      mask_size = length(mask)

      if (mask_size == 0) {
        #Update progress
        progress = progress + 1
        utils::setTxtProgressBar(pb, progress)
        next
      }

      beam_shot_n = level2a@h5[[beam_id]][["shot_number"]]$dims
      dt_dim = level2a@h5[[dt]]$dims

      if (length(dt_dim)[1] == 1) {
        if (dt_dim == 1) {
          hdf5r::createDataSet(newFile,dt,level2a@h5[[dt]][])
        } else if (dt_dim == beam_shot_n) {
          hdf5r::createDataSet(newFile,dt,level2a@h5[[dt]][mask])
        } else if ((dt_dim %% beam_shot_n) == 0) {
          n_waveforms = level2a@h5[[dt]]$dims / beam_shot_n
          v.seq = Vectorize(seq.default,vectorize.args = c("from"), SIMPLIFY=T)
          mask_init = mask*n_waveforms - (n_waveforms - 1)
          mask_waveform = matrix(v.seq(mask_init, len=n_waveforms), nrow=1)[1,]
          waveform=level2a@h5[[dt]][mask_waveform]
          hdf5r::createDataSet(newFile,dt,waveform)
        }
      } else if (length(dt_dim) == 2 && dt_dim[1] == beam_shot_n) {
        hdf5r::createDataSet(newFile,dt,level2a@h5[[dt]][mask,][])
      } else if (dt_dim[2] == beam_shot_n) {
        newFile$create_dataset(dt,level2a@h5[[dt]][1:dt_dim[1],mask])
      }
      else {
        stop(paste0("Don't know how to handle the dataset: ", dt, "\nContact the maintainer of the package!"))
      }

      #Update progress
      progress = progress + 1
      utils::setTxtProgressBar(pb, progress)
    }
  }

  hdf5r::h5flush(newFile)
  result = new("gedi.level2a", h5 = newFile)
  return (result)
}

#'Clip GEDI Level2A data by Coordinates
#'
#'@description This function clips GEDI Level2A data within given bounding coordinates
#'
#'
#'@param level2a A GEDI Level2A object (output of \code{\link[rGEDI:readLevel2A]{readLevel2A}} function). A S4 object of class "gedi.level2a".
#'@param xleft Numeric. West longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param xright Numeric. East longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param ybottom Numeric. South latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'@param ytopNumeric. North latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'@param output Optional character path where to save the new hdf5file. The default stores a temporary file only.
#'
#'@return An S4 object of class "gedi.level2a".
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi02_av001/
#'
#'@examples
#'# specify the path to GEDI level2A data
#'level2apath <- system.file("extdata", "GEDIexample_level02A.h5", package="rGEDI")
#'
#'# Reading GEDI level2A data
#'level2a<-readLevel2A(level1bpath)
#'
#'# Bounding rectangle coordinates
#'xleft = -44.15036
#'xright = -44.10066
#'ybottom = -13.75831
#'ytop = -13.71244
#'
#'# clip by extent boundary box
#'level2a_clip <- clipLevel2A(level1a,xleft,xright,ybottom,ytop)
#'
#'@export
clipLevel2A = function(level2a, xleft, xright, ybottom, ytop, output=""){
  if (output == "") {
    output = tempfile(fileext = ".h5")
  }
  output = fs::path_ext_set(output, "h5")

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
  newFile@h5$close_all()
  result = readLevel2A(output)

  return (result)
}

#'Clip GEDI Level2A data by geometry
#'
#'@description This function clips GEDI Level2A data within given geometry
#'
#'
#'@param level2a A GEDI Level2A object (output of \code{\link[rGEDI:readLevel2A]{readLevel2A}} function). A S4 object of class "gedi.level2a".
#'@param polygon_spdf Polygon. An object of class \code{\link[sp]{SpatialPolygonsDataFrame-class}},
#'which can be loaded as an ESRI shapefile using \code{\link[rgdal:readOGR]{readOGR}} function in the \emph{rgdal} package.
#'@param output optional character path where to save the new h5file. Default "" (temporary file).
#'
#'@return Returns a S4 object of class "gedi.level2a".
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi02_av001/
#'
#'@examples
#'
#'# specify the path to GEDI level2A data
#'level2apath <- system.file("extdata", "GEDIexample_level02A.h5", package="rGEDI")
#'
#'# Reading GEDI level2A data
#'level2a<-readLevel2A(level1bpath)
#'
#'# specify the path to shapefile
#'polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package="rGEDI")
#'
#'# Reading shapefile as SpatialPolygonsDataFrame object
#'library(rgdal)
#'polygon_spdf<-readOGR(polygons_filepath)
#'
#'level2a_clip <- clipLevel2AGeometry(level2a, polygon_spdf)
#'
#'@export
clipLevel2AGeometry = function(level2a, polygon_spdf, split_by = NULL, output="") {

  if (output == "") {
    output = tempfile(fileext = ".h5")
  }
  output = fs::path_ext_set(output, "h5")

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
      pts = suppressPackageStartupMessages(raster::intersect(points, polygon_spdf))

      mask_name = names(masks2)[i]
      if (is.null(split_by)) {
        polygon_masks[[""]][[beam]][[mask_name]] = pts@data[,1]
      } else {
        for (pol_id in as.character(unique(pts@data[split_by])[,1])) {
          polygon_masks[[pol_id]][[beam]][[mask_name]] = pts[(pts@data[split_by] == pol_id)[,1],]@data[,1]
        }
      }

      progress = progress + 1
      utils::setTxtProgressBar(pb, progress)
    }
  }
  close(pb)

  message("Writing new HDF5 files...")

  results = list()
  output="level2a_res.h5"
  i = 0
  len_masks = length(polygon_masks)
  for (pol_id in names(polygon_masks)) {
    i = i + 1
    message(gettextf("Writing %s='%s': %d of %d", split_by, pol_id, i, len_masks))
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
      h5_dt = level1b@h5[[dt]]
      dtype = h5_dt$get_type()
      if (is.na(all(h5_dt$chunk_dims))) {
        chunkdims = NULL
      } else {
        chunkdims = h5_dt$chunk_dims
      }

      if (length(dt_dim)[1] == 1) {
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
        hdf5r::createDataSet(newFile,dt,h5_dt[mask,][])
      } else if (dt_dim[2] == beam_shot_n) {
        newFile$create_dataset(dt,h5_dt[1:dt_dim[1],mask])
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

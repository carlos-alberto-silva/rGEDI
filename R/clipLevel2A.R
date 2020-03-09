#'Clip GEDI Level2A data by Coordinates
#'
#'@description This function clips GEDI Level2A data within given bounding coordinates
#'
#'@usage clipLevel2A(level2a, xmin, xmax, ymin, ymax, output)
#'
#'@param level2a A GEDI Level2A object (output of \code{\link[rGEDI:readLevel2A]{readLevel2A}} function). A S4 object of class "gedi.level2a".
#'@param xmin Numeric. West longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param xmax Numeric. East longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param ymin Numeric. South latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'@param ymax Numeric. North latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'@param output Optional character path where to save the new hdf5file. The default stores a temporary file only.
#'
#'@return Returns a list of S4 object of class "gedi.level2a".
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi02_av001/
#'
#'@examples
#'# specify the path to GEDI level2A data (zip file)
#'level2A_fp_zip <- system.file("extdata",
#'                   "GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level2A data
#'level2Apath <- unzip(level2A_fp_zip,exdir = dirname(level2A_fp_zip))
#'
#'# Reading GEDI level2A data (h5 file)
#'level2a<-readLevel2A(level2Apath=level2Apath)
#'
#'# Bounding rectangle coordinates
#'xmin = -44.15036
#'xmax = -44.10066
#'ymin = -13.75831
#'ymax = -13.71244
#'
#'# Spepecifing output file and path
#'output<-file.path(getwd(),"GEDI02_A_2019108080338_O01964_T05337_02_001_01_clip.h5")
#'
#'# clip by extent boundary box
#'level2a_clip <- clipLevel2A(level2a,xmin,xmax,ymin,ymax,output)
#'
#'close(level2a)
#'close(level2a_clip)
#'@export
clipLevel2A = function(level2a, xmin, xmax, ymin, ymax, output=""){
  if (output == "") {
    output = tempfile(fileext = ".h5")
  }
  output = fs::path_ext_set(output, "h5")

  # Get all spatial data as a list of dataframes with spatial information
  spData = getSpatialData2A(level2a)

  masks = clipSpDataByExtentLevel2A(spData, xmin, xmax, ymin, ymax)

  newFile = clipByMask2A(level2a,
                         masks,
                         output)
  output = newFile@h5$filename
  close(newFile)
  result = readLevel2A(output)

  return (result)
}

#'Clip GEDI Level2A data by geometry
#'
#'@description This function clips GEDI Level2A data within given geometry
#'
#'@usage clipLevel2AGeometry(level2a, polygon_spdf, output="", split_by=NULL)
#'
#'@param level2a A GEDI Level2A object (output of \code{\link[rGEDI:readLevel2A]{readLevel2A}} function). A S4 object of class "gedi.level2a".
#'@param polygon_spdf Polygon. An object of class \code{\link[sp]{SpatialPolygonsDataFrame-class}},
#'which can be loaded as an ESRI shapefile using \code{\link[raster:shapefile]{raster::shapefile()}} function in the \emph{raster} package.
#'@param output optional character path where to save the new h5file. Default "" (temporary file).
#'@param split_by Polygon id. If defined, GEDI data will be clipped by each polygon using attribute specified by \code{split_by} from attribute table.
#'
#'@return Returns a list of S4 object of class "gedi.level2a".
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi02_av001/
#'
#'@examples
#'# specify the path to GEDI level2A data (zip file)
#'level2A_fp_zip <- system.file("extdata",
#'                   "GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level2A data
#'level2Apath <- unzip(level2A_fp_zip,exdir = dirname(level2A_fp_zip))
#'
#'# Reading GEDI level2A data (h5 file)
#'level2a<-readLevel2A(level2Apath=level2Apath)
#'
#'# specify the path to shapefile
#'polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package="rGEDI")
#'
#'# Reading shapefile as SpatialPolygonsDataFrame object
#'library(raster)
#'polygon_spdf<-shapefile(polygon_filepath)
#'
#'# Spepecifing output file and path
#'output<-file.path(getwd(),"GEDI02_A_2019108080338_O01964_T05337_02_001_01_clip")
#'
#'# clip by extent boundary box
#'level2a_clip <- clipLevel2AGeometry(level2a, polygon_spdf = polygon_spdf,
#'                                    output=output,
#'                                    split_by="id")
#'close(level2a)
#'lapply(level2a_clip, close)
#'@export
clipLevel2AGeometry = function(level2a, polygon_spdf, output="", split_by = NULL) {
  output = checkOutput(output)
  stopifnotMessage(
    "split_by not in polygon_spdf"=is.null(split_by) || split_by %in% colnames(polygon_spdf@data)
  )

  spData = getSpatialData2A(level2a)

  xmin = polygon_spdf@bbox[1,1]
  xmax = polygon_spdf@bbox[1,2]
  ymin = polygon_spdf@bbox[2,1]
  ymax = polygon_spdf@bbox[2,2]

  masks = clipSpDataByExtentLevel2A(spData, xmin, xmax, ymin, ymax)

  polygon_masks = getPolygonMaskLevel2A(spData, masks, polygon_spdf, split_by)

  results = clipByMasks(level2a, polygon_masks, output, split_by, clipByMask2A)

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
  if(length(masks) == 0) {
    message("No intersection found!")
    newFile$close_all()
    newFile = hdf5r::H5File$new(output, mode="r")
    result = new("gedi.level2a", h5 = newFile)
    return(result)
  }

  for (attr in hdf5r::list.attributes(level2a@h5)) {
    hdf5r::h5attr(newFile, attr) = hdf5r::h5attr(level2a@h5, attr)
  }


  all_groups = hdf5r::list.groups(level2a@h5)


  # Check if the beam has any intersecting area
  beams_with_value = lapply(lapply(masks, function(x) sapply(x, length)), sum)>0
  beams_with_value = names(which(beams_with_value))
  beams_with_value = c(beams_with_value, "METADATA")
  which_groups = gsub("([^/]*).*","\\1",all_groups) %in% beams_with_value
  groups_with_value = all_groups[which_groups]

  # Setup progress bar
  all_datasets = hdf5r::list.datasets(level2a@h5)
  which_datasets = gsub("([^/]*).*","\\1",all_datasets) %in% beams_with_value
  datasets_with_value = all_datasets[which_datasets]

  total = length(datasets_with_value)
  pb = utils::txtProgressBar(min = 0, max = total, style = 3)
  progress = 0

  for (group in groups_with_value) {
    beam_id = strsplit(group, "/")[[1]][1]

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
        beam_shot_n = level2a@h5[[beam_id]][["shot_number"]]$dims
      } else {
        mask = masks[[beam_id]][[alg_id]]
        beam_shot_n = level2a@h5[[sprintf("%s/geolocation/elev_highestreturn_a%d", beam_id, alg_id)]]$dims
      }
      mask_size = length(mask)

      if (mask_size == 0) {
        #Update progress
        progress = progress + 1
        utils::setTxtProgressBar(pb, progress)
        next
      }


      h5_dt = level2a@h5[[dt]]
      dt_dim = h5_dt$dims
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
        if (length(mask) == 1) {
          chunkdims = chunkdims[[2]]
        }
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
  newFile$close_all()
  newFile =  hdf5r::H5File$new(output, mode="r")
  result = new("gedi.level2a", h5 = newFile)
  close(pb)
  return (result)
}

#'Clip GEDI Level1B data by Coordinates
#'
#'@description This function clips GEDI Level1B data within given bounding coordinates
#'
#'@usage clipLevel1B(level1b, xmin, xmax, ymin, ymax, output)
#'
#'@param level1b A GEDI Level1B object (output of \code{\link[rGEDI:readLevel1B]{readLevel1B}} function). A S4 object of class "gedi.level1b".
#'@param xmin Numeric. West longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param xmax Numeric. East longitude (x) coordinate of bounding rectangle, in decimal degrees.
#'@param ymin Numeric. South latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'@param ymax Numeric. North latitude (y) coordinate of bounding rectangle, in decimal degrees.
#'@param output Optional character path where to save the new hdf5file. The default stores a temporary file only.
#'
#'@return An S4 object of class "gedi.level1b".
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi01_bv001/
#'
#'@examples
#'\dontrun{
#'# specify the path to download GEDI example dataset
#'outdir<-getwd()
#'
#'# downloading GEDI example dataset (zip file)
#'download.file(
#'              paste0(
#'                     "https://github.com/carlos-alberto-silva/rGEDI/",
#'                     "releases/download/examples/examples.zip"
#'              ),
#'              destfile=paste0(outdir,"/examples.zip"))
#'
#'# unzip the file
#'unzip(paste0(outdir,"\\examples.zip"))
#'
#'# specify the path to GEDI level1B data
#'level1bpath = paste0(outdir,"\\GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.h5"))
#'
#'# Reading GEDI level1B file
#'level1b<-readLevel1b(gedilevel1b)
#'
#'# Bounding rectangle coordinates
#'xmin = -44.15036
#'xmax = -44.10066
#'ymin = -13.75831
#'ymax = -13.71244
#'
#'# clip by extent boundary box
#'level1b_clip <- clipLevel1B(level1b,xmin, xmax, ymin, ymax)
#'
#'}
#'@import hdf5r fs
#'@export
clipLevel1B = function(level1b, xmin, xmax, ymin, ymax, output=""){
  output = checkOutput(output)

  # Get all spatial data as a list of dataframes with spatial information
  spData = getSpatialData1B(level1b)

  masks = clipSpDataByExtentLevelB(spData, xmin, xmax, ymin, ymax)

  newFile = clipByMask1B(level1b,
                         masks,
                         output)
  output = newFile@h5$filename
  newFile@h5$close_all()
  result = readLevel1B(output)

  return (result)
}

#'Clip GEDI Level1B data by geometry
#'
#'@description This function clips GEDI Level1B data within a bounding geometry
#'
#'
#'@param level1b A GEDI Level1B object (output of \code{\link[rGEDI:readLevel1B]{readLevel1B}} function). A S4 object of class "gedi.level1b".
#'@param polygon_spdf Polygon. An object of class \code{\link[sp]{SpatialPolygonsDataFrame-class}},
#'which can be loaded as an ESRI shapefile using \code{\link[rgdal:readOGR]{readOGR}} function in the \emph{rgdal} package.
#'@param split_by Polygon id. If defined, GEDI data will be clipped by each polygon using the polygon id from table of attribute defined by the user
#'@param output Optional character path where to save the new hdf5file. The default stores a temporary file only.
#'
#'@return An S4 object of class "gedi.level1b".
#'
#'@examples
#'\dontrun{
#'# specify the path to download GEDI example dataset
#'outdir<-getwd()
#'
#'# downloading GEDI example dataset (zip file)
#'download.file(
#'              paste0(
#'                     "https://github.com/carlos-alberto-silva/rGEDI/",
#'                     "releases/download/examples/examples.zip"
#'              ),
#'              destfile=paste0(outdir,"/examples.zip"))
#'
#'# unzip the file
#'unzip(paste0(outdir,"\\examples.zip"))
#'
#'# specify the path to GEDI level1B data
#'level1bpath = paste0(outdir,"\\GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.h5")
#'
#'# Reading GEDI level1B file
#'level1b<-readLevel1b(gedilevel1b)
#'
#'# specify the path to shapefile
#'polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package="rGEDI")
#'
#'# Reading shapefile as SpatialPolygonsDataFrame object
#'library(rgdal)
#'polygon_spdf<-readOGR(polygons_filepath)
#'
#'# clip by extent boundary box
#'level1b_clip <- clipLevel1BGeometry(level1b, polygon_spdf = polygon_spdf)
#'}
#'@import hdf5r
#'@export
clipLevel1BGeometry = function(level1b, polygon_spdf, output="", split_by=NULL) {
  output = checkOutput(output)

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

    hdf5r::createGroup(newFile, group)
    createAttributesWithinGroup(level1b@h5, newFile, group)

    for (dt in hdf5r::list.datasets(level1b@h5[[group]], recursive = FALSE, full.names = T)) {
      beam_shot_n = level1b@h5[[beam_id]][["shot_number"]]$dims
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
        } else if ((dt_dim %% beam_shot_n) == 0) {
          n_waveforms = level1b@h5[[dt]]$dims / beam_shot_n
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
        hdf5r::createDataSet(newFile,dt,level1b@h5[[dt]][mask,], dtype=dtype, chunk_dim=chunkdims)
      } else {
        stop(paste0("Don't know how to treat dataset: ", dt, "\nContact the maintainer of the package!"))
      }

      #Update progress
      progress = progress + 1
      utils::setTxtProgressBar(pb, progress)
    }
  }

  level1b@h5 = newFile
  close(pb)
  #spatial = level1B2dt(level1b)
  return (level1b)
}

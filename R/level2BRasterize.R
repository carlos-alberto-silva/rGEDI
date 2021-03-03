#' Aggregate selected metrics into raster tif cells
#'
#' @description
#' This function will read multiple Level2B H5 files and aggregate into multiple rasters: count, and 1st, 2nd, 3rd and 4th moments (count, m1, m2, m3 and m4) for each metric selected, from which we can calculate statistics such as Mean, SD, Skewness and Kurtosis.
#'
#' @param l2bDir CharacterVector. The directory paths where the H5 GEDI files are stored;
#' @param metrics CharacterVector. A vector of metrics available from Level2B, as in the getLevel2BVPM documentation
#' @param out_root Character. The root name for the raster output files, the pattern is {out_root}_{metric}_{count/m1/m2/m3/m4}.tif. This should include the full path for the file.
#' @param ul_lat Numeric. Upper left latitude for the bounding box
#' @param ul_lon Numeric. Upper left longitude for the bounding box
#' @param lr_lat Numeric. Lower right latitude for the bounding box
#' @param lr_lon Numeric. Lower right longitude for the bounding box
#' @param res NumericVector. Resolution lon lat for the output raster in coordinates decimal degrees
#' @param creation_options CharacterVector. The GDAL creation options for the tif file. Default c("COMPRESS=PACKBITS", "BIGTIFF=IF_SAFER", "TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512") will create BIGTIFF if needed, with PACKBITS compression and tiled by 512x512 pixels.
#'
#' @details 
#' This function will create five different aggregate statistics (count, m1, m2, m3 and m4). m1 to m4 are the central moments. One can calculate mean, standard deviation, skewness and kurtosis with the following formulas according to Terriberry (2007) and \insertCite{Joanes1998;textual}{rGEDI}:
#'
#' \deqn{ \bar{x} = m_1 }{mean = m1}
#' 
#' \deqn{ s = \sqrt{\frac{m_2}{count - 1}} }{sd = sqrt(m2/(count - 1))}
#' 
#' \deqn{ g_1 = \frac{\sqrt{count} \cdot m_3}{m_2^{1.5}} }{ g1 = (sqrt(count) * m3) / (m2^1.5)}
#' 
#' \deqn{ g_2 = \frac{count \cdot m_4}{m_2^2} - 3 }{g2 = (count * m4) / (m2 * m2) - 3.0}
#' 
#' \deqn{ skewness = \frac{\sqrt{count(count - 1)}}{n-2} g_1 }{skewness = sqrt((count * (count - 1))) * g1 / (count - 2)}
#' 
#' \deqn{ kurtosis = \frac{count - 1}{(count - 2)(count - 3)}[(count + 1)g_2 + 6] }{kurtosis = ((count - 1) / ((count - 2) * (count - 3))) * ((count + 1) * g2 + 6)}
#' 
#' @references
#' \insertAllCited{}
#' 
#' Terriberry, Timothy B. (2007), Computing Higher-Order Moments Online, archived from the original on 23 April 2014, retrieved 5 May 2008
#' 
#' @return Nothing. It outputs multiple raster tif files to the out_root specified path.
#'
#' @examples
#' # Specifying the path to GEDI level2B data (zip file)
#' outdir = tempdir()
#' level2B_fp_zip <- system.file("extdata",
#'                    "GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'                    package="rGEDI")
#'
#' # Unzipping GEDI level2A data
#' level2Bpath <- unzip(level2B_fp_zip,exdir = outdir)
#'
#' # Reading GEDI level2B data (h5 file)
#' level2b<-readLevel2B(level2Bpath=level2Bpath)
#'
#' ul_lat = -13.72016
#' ul_lon = -44.14000
#' lr_lat = -13.74998
#' lr_lon = -44.11009
#'
#' res = 100 # meters
#' lat_to_met_factor = 1 / 110540
#' lon_to_met_factor = 1 / 111320
#' xres = lon_to_met_factor * res
#' yres = lat_to_met_factor * res
#'
#' rast = level2bRasterizeStats(
#'   outdir,
#'   metrics = c("rh100"),
#'   out_root = file.path(outdir, "output"),
#'   ul_lat = -13.72016,
#'   ul_lon = -44.14000,
#'   lr_lat = -13.74998,
#'   lr_lon = -44.11009,
#'   res = c(xres, yres)
#'   )
#'
#' close(level2b)
#'
#' @import e1071
#' @import data.table
#' @export
level2bRasterizeStats = function(l2bDir,
    metrics, out_root, ul_lat, ul_lon, lr_lat, lr_lon,
    res, creation_options = c("COMPRESS=PACKBITS",
        "BIGTIFF=IF_SAFER",
        "TILED=YES",
        "BLOCKXSIZE=512",
        "BLOCKYSIZE=512"
        )
    ) {
        x_blocks =
        y_blocks =
        l2b_quality_flag =
        cols_without_quality =
        x_ind =
        longitude_bin0 =
        y_ind =
        latitude_bin0 =
        inds = NULL
  projstring = 'GEOGCS["WGS 84",
    DATUM["WGS_1984",
        SPHEROID["WGS 84",6378137,298.257223563,
            AUTHORITY["EPSG","7030"]],
        AUTHORITY["EPSG","6326"]],
    PRIMEM["Greenwich",0,
        AUTHORITY["EPSG","8901"]],
    UNIT["degree",0.01745329251994328,
        AUTHORITY["EPSG","9122"]],
    AUTHORITY["EPSG","4326"]]'
  l2b_list = sapply(l2bDir, function(search_path) list.files(search_path, "GEDI02_B.*h5"))
  total_files = length(l2b_list)

  lat_min = lr_lat
  lon_min = ul_lon
  lat_max = ul_lat
  lon_max = lr_lon
  xres = res[1]
  yres = res[2]


  cols.coord = c("latitude_bin0", "longitude_bin0", "l2b_quality_flag")

  metricCounter = 0
  nMetrics = length(metrics)
  # metric = metrics[1]
  for (metric in metrics) {
    metricCounter = metricCounter + 1
    message(sprintf("Metric %s (%d/%d)", metric, metricCounter, nMetrics), appendLF = T)
    cols = c(cols.coord, metric)

    count_path = sprintf("%s_%s_%s.tif", out_root, metric, "count")
    m1_path = sprintf("%s_%s_%s.tif", out_root, metric, "m1")
    m2_path = sprintf("%s_%s_%s.tif", out_root, metric, "m2")
    m3_path = sprintf("%s_%s_%s.tif", out_root, metric, "m3")
    m4_path = sprintf("%s_%s_%s.tif", out_root, metric, "m4")

    count_rast = GDALDataset$new(
    raster_path = count_path,
    nbands = 1,
    datatype = GDALDataType$GDT_Int32,
    projstring = projstring,
    lr_lat = lr_lat,
    ul_lat = ul_lat,
    ul_lon = ul_lon,
    lr_lon = lr_lon,
    res = c(xres, - yres),
    nodata = 0,
    co = creation_options)

    xsize = count_rast$GetRasterXSize()
    ysize = count_rast$GetRasterYSize()

    m1_rast = GDALDataset$new(
    raster_path = m1_path,
    nbands = 1,
    datatype = GDALDataType$GDT_Float32,
    projstring = projstring,
    lr_lat = lr_lat,
    ul_lat = ul_lat,
    ul_lon = ul_lon,
    lr_lon = lr_lon,
    res = c(xres, - yres),
    nodata = -9999,
    co = creation_options)

    m2_rast = GDALDataset$new(
    raster_path = m2_path,
    nbands = 1,
    datatype = GDALDataType$GDT_Float32,
    projstring = projstring,
    lr_lat = lr_lat,
    ul_lat = ul_lat,
    ul_lon = ul_lon,
    lr_lon = lr_lon,
    res = c(xres, - yres),
    nodata = -9999,
    co = creation_options)

    m3_rast = GDALDataset$new(
    raster_path = m3_path,
    nbands = 1,
    datatype = GDALDataType$GDT_Float32,
    projstring = projstring,
    lr_lat = lr_lat,
    ul_lat = ul_lat,
    ul_lon = ul_lon,
    lr_lon = lr_lon,
    res = c(xres, - yres),
    nodata = -9999,
    co = creation_options)

    m4_rast = GDALDataset$new(
    raster_path = m4_path,
    nbands = 1,
    datatype = GDALDataType$GDT_Float32,
    projstring = projstring,
    lr_lat = lr_lat,
    ul_lat = ul_lat,
    ul_lon = ul_lon,
    lr_lon = lr_lon,
    res = c(xres, - yres),
    nodata = -9999,
    co = creation_options)


    n_band = count_rast[[1]]
    m1_band = m1_rast[[1]]
    m2_band = m2_rast[[1]]
    m3_band = m3_rast[[1]]
    m4_band = m4_rast[[1]]

    block_x_size = n_band$GetBlockXSize()
    block_y_size = n_band$GetBlockYSize()

    file_index = 0
    # l2b_path = l2b_list[1]
    for (l2b_path in l2b_list) {
      file_index = file_index + 1
      message(sprintf("Reading file %s (%d/%d)", l2b_path, file_index, total_files), appendLF = T)
      l2b = readLevel2B(file.path(l2bDir, l2b_path))

      vals = getLevel2BVPM(l2b, cols = cols)
      cols_without_quality = c(setdiff(cols, "l2b_quality_flag"))
      vals = vals[l2b_quality_flag == 1, ..cols_without_quality]


      vals[, x_ind := as.integer(vals[, floor((longitude_bin0 - lon_min) / xres)])]
      vals[, y_ind := as.integer(vals[, 1+floor((latitude_bin0 - lat_min) / yres)])]

      blocks = vals[, lapply(.SD, function(x) as.integer(floor(x / block_x_size))), .SDcols = c("x_ind", "y_ind")]
      colnames(blocks) = c("x_blocks", "y_blocks")
      vals = cbind(vals, blocks)
      df_unique = unique(blocks)
      df_unique = df_unique[x_blocks >= 0 & x_blocks <= floor(xsize / block_x_size) & y_blocks >= 0 & y_blocks <= floor(ysize / block_y_size)]
      total_rows = nrow(df_unique)
      # ii = 1
      for (ii in 1:total_rows) {
        message(sprintf("\rProcessing blocks...%2f", (100.0*ii)/total_rows), appendLF = F)
        row = df_unique[ii,]
        x_block = row$x_blocks
        y_block = row$y_blocks
        this_vals = vals[x_blocks == x_block & y_blocks == y_block]
        agg1 = data.table::data.table(
        n = n_band$ReadBlock(x_block, y_block),
        M1 = m1_band$ReadBlock(x_block, y_block),
        M2 = m2_band$ReadBlock(x_block, y_block),
        M3 = m3_band$ReadBlock(x_block, y_block),
        M4 = m4_band$ReadBlock(x_block, y_block)
      )

        this_vals[, x_ind := x_ind - x_block * block_x_size]
        this_vals[, y_ind := y_ind - y_block * block_y_size]
        this_vals[, inds := x_ind + y_ind * block_x_size]

        agg2 = this_vals[, list(
        n = length(eval(as.name(metric))),
        M1 = mean(eval(as.name(metric))),
        M2 = e1071::moment(eval(as.name(metric)), order = 2, center = TRUE) * length(eval(as.name(metric))),
        M3 = e1071::moment(eval(as.name(metric)), order = 3, center = TRUE) * length(eval(as.name(metric))),
        M4 = e1071::moment(eval(as.name(metric)), order = 4, center = TRUE) * length(eval(as.name(metric)))
      ), by = list(inds)]
        agg1[agg2$inds] = combine_stats(agg1[agg2$inds], agg2)

        n_band$WriteBlock(x_block, y_block, agg1$n)
        m1_band$WriteBlock(x_block, y_block, agg1$M1)
        m2_band$WriteBlock(x_block, y_block, agg1$M2)
        m3_band$WriteBlock(x_block, y_block, agg1$M3)
        m4_band$WriteBlock(x_block, y_block, agg1$M4)
      }
      close(l2b)
    }
  }
  gc()
}

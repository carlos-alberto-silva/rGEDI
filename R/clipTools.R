#'@import sp
NULL

createAttributesWithinGroup = function(h5, newFile, group = "/") {
  for (attr in hdf5r::list.attributes(h5[[group]])) {
    hdf5r::h5attr(newFile[[group]], attr) = hdf5r::h5attr(h5[[group]], attr)
  }
}

clipSpDataByExtentLevel2A = function(spData, xmin, xmax, ymin, ymax) {
  masks = lapply(spData, function(x) {
    masks2 = lapply(x, function(y) {
      mask = y$longitude_lowest >= xmin &
        y$longitude_lowest <= xmax &
        y$latitude_highest >= ymin &
        y$latitude_highest <= ymax &
        y$longitude_highest >= xmin &
        y$longitude_highest <= xmax &
        y$latitude_lowest >= ymin &
        y$latitude_lowest <= ymax

      mask[!stats::complete.cases(mask)] = FALSE
      return((seq_along(y$longitude_lowest))[mask])
    })
    return (masks2)
  })
  if(all(sapply(masks, function(x) sum(sapply(x, length)))==0)){
    stop("The clipping ROI does not intersect with the data!")
  }
  return(masks)
}


getPolygonMaskLevel2A = function(spData, masks, polygon, split_by) {
  message("Intersecting with polygon...")
  pb = utils::txtProgressBar(min = 0, max = length(masks), style = 3)
  progress = 0
  polygon_masks = list()

  if (is.null(split_by)) {
    masknames = ""
  } else {
    masknames = unique(paste0(polygon[[split_by]]))
  }

  for (m in masknames) {
    polygon_masks[[m]] = list()
  }

  for (beam in names(masks)) {
    masks2 = masks[[beam]]

    for (i in seq_along(masks2)) {
      mask = masks2[[i]]
      if (length(mask) == 0) next

      spDataMasked = spData[[beam]][[i]][mask, ]
      points = sf::st_as_sf(
        spDataMasked,
        coords = c("longitude_highest", "latitude_highest"),
        crs = sf::st_crs(polygon)
      )

      pts = suppressPackageStartupMessages(raster::intersect(points, polygon))


      mask_name = names(masks2)[i]
      if (ncol(pts@data) == 2) {
        split_by2 = 2
      } else {
        split_by2 = split_by
      }
      if (is.null(split_by)) {
        polygon_masks[[1]][[beam]][[mask_name]] = pts[, 1]
      } else {
        for (pol_id in as.character(unique(pts[split_by2])[, 1])) {
          polygon_masks[[pol_id]][[beam]][[mask_name]] = pts[(pts[split_by2] == pol_id)[, 1], ][, 1]
        }
      }

      progress = progress + 1
      utils::setTxtProgressBar(pb, progress)
    }
  }
  close(pb)

  if (all(
    sapply(
      polygon_masks, function(x) {
        res=sapply(x, function(y) sum(sapply(y, length)))
        ifelse(length(res)>0, sum(res), 0)
      }
    ) == 0))
  {
    stop("The clipping polygon does not intersect with the data!")
  }
  return(polygon_masks)
}


clipSpDataByExtentLevelB = function(spData, xmin, xmax, ymin, ymax) {
  masks = lapply(spData, function(x) {
    mask = x$longitude_bin0 >= xmin &
      x$longitude_bin0 <= xmax &
      x$latitude_bin0 >= ymin &
      x$latitude_bin0 <= ymax &
      x$longitude_lastbin >= xmin &
      x$longitude_lastbin <= xmax &
      x$latitude_lastbin >= ymin &
      x$latitude_lastbin <= ymax

    mask[!stats::complete.cases(mask)] = FALSE
    return((seq_along(x$longitude_bin0))[mask])
  })
  if (all(sapply(masks, length)==0)) {
    stop("The clipping ROI does not intersect with the data!")
  }
  return (masks)
}

checkOutput = function(output) {
  if (output == "") {
    output = tempfile(fileext = ".h5")
  }
  output = fs::path_ext_set(output, "h5")
  return(output)
}


getPolygonMaskLevelB = function(spData, masks, polygon, split_by) {
  message("Intersecting with polygons...")
  pb = utils::txtProgressBar(min = 0, max = length(masks), style = 3)
  progress = 0
  polygon_masks = list()

  if (is.null(split_by)) {
    masknames = ""
  } else {
    masknames = unique(paste0(polygon[[split_by]]))
  }
  for (m in masknames) {
    polygon_masks[[m]] = list()
  }

  for (beam in names(masks)) {
    mask = masks[[beam]]

    if (length(mask) == 0) next

    spDataMasked = spData[[beam]][mask, ]
    points = sf::st_as_sf(
      spDataMasked,
      coords = c("longitude_bin0", "latitude_bin0"),
      crs = sf::st_crs(polygon)
    )

    pts = suppressPackageStartupMessages(terra::intersect(points, polygon))
    if (ncol(pts) == 2) {
      split_by2 = 2
    } else {
      split_by2 = split_by
    }
    if (is.null(split_by)) {
      polygon_masks[[1]][[beam]] = pts[,1]
    } else {
      for (pol_id in unique(as.character(paste0(pts[[split_by2]])))) {

        polygon_masks[[pol_id]][[beam]] = pts[pts[[split_by2]] == pol_id,][,1]
      }
    }

    progress = progress + 1
    utils::setTxtProgressBar(pb, progress)
  }
  close(pb)

  if (all(sapply(polygon_masks, length)==0)) {
    stop("The polygon does not intersect with the data!")
  }
  return(polygon_masks)
}

clipByMasks = function(h5file, polygon_masks, output, split_by, clipFun) {
  message("Writing new HDF5 files...")
  results = list()
  i = 0
  len_masks = length(polygon_masks)
  for (pol_idx in seq_along(polygon_masks)) {
    pol_id = names(polygon_masks)[pol_idx]
    i = i + 1
    message(sprintf("Writing %s='%s': %d of %d", split_by, pol_id, i, len_masks))
    output2 = gsub("\\.h5$", paste0("_", pol_id,".h5"), output)
    results[[pol_id]] = clipFun(h5file,
                                masks = polygon_masks[[pol_idx]],
                                output = output2)
  }

  return(results)
}


checkClipExtentInputs = function(obj, className, xmin, xmax, ymin, ymax) {
  criterias = list()
  criterias[paste0("Object is not from class", className)] = class(obj) == className
  criterias = c(criterias, list(
    "xmin is not numeric" = class(xmin) == "numeric",
    "xmax is not numeric" = class(xmax) == "numeric",
    "ymin is not numeric" = class(ymin) == "numeric",
    "ymax is not numeric" = class(ymax) == "numeric"
  ))
  do.call(stopifnotMessage, criterias)
}

checkClipGeoInputs = function(obj, className, polygon, split_by) {
  criterias = list()
  criterias[paste0("Object is not from class", className)] = class(obj) == className
  criterias = c(criterias, list(
    "is not a valid sf object" = all(c('sf', 'data.frame') %in% class(polygon)),
    "polygon_spdf is not a SpatialPolygonsDataFrame" = any(c('POLYGON', 'MULTIPOLYGON') %in% sf::st_geometry_type(polygon)),
    "split_by is not a valid attribute of polygon_spdf" = is.null(split_by) || split_by %in% colnames(polygon)
  ))
  do.call(stopifnotMessage, criterias)
}

#!/usr/bin/env Rscript

# Build the complete rGEDI public-function map from NAMESPACE.
#
# The figure is intentionally generated with base R's grid graphics so it can
# be rebuilt on Windows, macOS, Linux, CI, and CRAN machines without optional
# plotting packages. The script stops if an exported function is missing from
# the grouping below or appears more than once.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1L]]) else "tools/build-function-map.R"
root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
namespace_file <- file.path(root, "NAMESPACE")
output_dir <- file.path(root, "output", "pdf")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

main_groups <- list(
  "1. Configure, find and access" = c(
    "rGEDI_configure", "gedifinder", "gediDownload", "earthdata_login",
    "ee_initialize"
  ),
  "2. Read products" = c(
    "readLevel1B", "readLevel2A", "readLevel2B", "readLevel3",
    "readLevel4A", "readLevel4B"
  ),
  "3. Get measurements" = c(
    "getLevel1BGeo", "getLevel1BWF", "getLevel2AM",
    "getLevel2BPAIProfile", "getLevel2BPAVDProfile", "getLevel2BVPM",
    "getLevel3", "getLevel4A", "getLevel4B"
  ),
  "4. Spatial statistics and grids" = c(
    "rasterize", "polyStats", "to_vect", "gridStatsLevel2AM",
    "gridStatsLevel2BVPM", "gridStatsLevel4A", "polyStatsLevel2AM",
    "polyStatsLevel2BVPM", "polyStatsLevel4A"
  ),
  "5. Clip products" = c(
    "clip", "clipLevel1B", "clipLevel1BGeo", "clipLevel1BGeoGeometry",
    "clipLevel1BGeometry", "clipLevel2A", "clipLevel2AGeometry",
    "clipLevel2AM", "clipLevel2AMGeometry", "clipLevel2B",
    "clipLevel2BGeometry", "clipLevel2BPAIProfile",
    "clipLevel2BPAIProfileGeometry", "clipLevel2BPAVDProfile",
    "clipLevel2BPAVDProfileGeometry", "clipLevel2BVPM",
    "clipLevel2BVPMGeometry", "clipLevel3", "clipLevel4A",
    "clipLevel4AGeometry", "clipLevel4B"
  ),
  "6. Sample, model and map" = c(
    "sample", "fit_model", "varSel", "fit_metrics", "map_create",
    "map_view", "map_download"
  ),
  "7. Simulate and visualize" = c(
    "gediWFSimulator", "gediWFMetrics", "getTrack",
    "plot_gedi_orbit_animation", "plotWFMetrics", "plotPAIProfile",
    "plotPAVDProfile", "plotLevel3", "plotLevel4A", "plotLevel4B"
  )
)

auxiliary_groups <- list(
  "Sampling specifications" = c(
    "randomSampling", "spacedSampling", "gridSampling",
    "stratifiedSampling", "geomSampling", "rasterSampling"
  ),
  "Raster product operations" = c(
    "mosaicLevel3", "mosaicLevel4B", "aggregateLevel3",
    "aggregateLevel4B", "projectLevel3", "projectLevel4B",
    "resampleLevel3", "resampleLevel4B"
  ),
  "GEE catalog, display and export" = c(
    "gediEE", "search_datasets", "get_catalog_id", "addEEImage",
    "getTileUrl", "ee_image_to_drive", "ee_check_task_status"
  ),
  "Cloud and package configuration" = c(
    "earthdata_s3_credentials", "gedidb_available", "gedidb_disconnect",
    "gee_available"
  ),
  "Archive and batch processing" = c(
    "accumulateGEDIGranule", "accumulateGEDIGranules",
    "assembleGEDIArchiveYears", "readGEDIArchiveRegion"
  ),
  "GEE predictors and geometry" = c(
    "build_ee_forest", "extractEE", "extractGEDIAncillary",
    "ee_build_AlphaEarth_embedding_terrain_stack",
    "ee_build_hls_s1c_terrain_stack", "ee_rect_to_sf", "ext_to_ee",
    "vect_as_ee", "slope", "aspect", "glcmTexture"
  ),
  "Classes and storage helpers" = c(
    "gedi.fullwaveform", "gedi.level1b", "gedi.level2a", "gedi.level2b",
    "gedi.level4a", "predictGEDIH5"
  )
)

namespace <- readLines(namespace_file, warn = FALSE)
exports <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", namespace, value = TRUE))
grouped <- unname(unlist(c(main_groups, auxiliary_groups), use.names = FALSE))
duplicates <- sort(unique(grouped[duplicated(grouped)]))
missing <- sort(setdiff(exports, grouped))
extra <- sort(setdiff(grouped, exports))
if (length(duplicates) || length(missing) || length(extra)) {
  stop(
    "Function-map grouping mismatch.\n",
    "Duplicates: ", paste(duplicates, collapse = ", "), "\n",
    "Missing: ", paste(missing, collapse = ", "), "\n",
    "Extra: ", paste(extra, collapse = ", ")
  )
}

s3_lines <- grep("^S3method\\(", namespace, value = TRUE)
s3_methods <- vapply(s3_lines, function(x) {
  fields <- strsplit(sub("^S3method\\((.*)\\)$", "\\1", x), ",", fixed = TRUE)[[1L]]
  generic <- gsub('^"|"$', "", trimws(fields[[1L]]))
  class <- trimws(fields[[2L]])
  paste0(generic, ".", class)
}, character(1L), USE.NAMES = FALSE)

colors <- list(
  background = "#F7F5F0", ink = "#173730", green = "#0C463D",
  green2 = "#176B54", pale = "#DCECE6", line = "#A9C5BD",
  orange = "#DE6327", muted = "#55716B", white = "#FFFFFF"
)

draw_box <- function(x, y, width, height, fill, border = colors$line,
                     radius = 0.008, lwd = 1.4) {
  grid::grid.roundrect(
    x = x, y = y, width = width, height = height,
    just = c("left", "bottom"), r = grid::unit(radius, "npc"),
    gp = grid::gpar(fill = fill, col = border, lwd = lwd)
  )
}

draw_card <- function(title, functions, x, y, width, height,
                      primary = TRUE, columns = 1L, predict = FALSE) {
  fill <- if (primary) colors$white else colors$pale
  border <- if (primary) colors$green2 else colors$line
  draw_box(x, y, width, height, fill, border, lwd = if (primary) 1.7 else 1.1)
  header_h <- if (primary) 0.030 else 0
  if (primary) {
    grid::grid.roundrect(
      x = x, y = y + height - header_h, width = width, height = header_h,
      just = c("left", "bottom"), r = grid::unit(0.007, "npc"),
      gp = grid::gpar(fill = colors$green2, col = colors$green2)
    )
    grid::grid.rect(
      x = x, y = y + height - header_h, width = width, height = header_h / 2,
      just = c("left", "bottom"),
      gp = grid::gpar(fill = colors$green2, col = colors$green2)
    )
  }
  grid::grid.text(
    title, x = x + 0.006, y = y + height - if (primary) 0.015 else 0.012,
    just = c("left", "center"),
    gp = grid::gpar(
      col = if (primary) colors$white else colors$green,
      fontfamily = "sans", fontface = "bold", fontsize = if (primary) 9 else 8.5
    )
  )
  n <- length(functions)
  rows <- ceiling(n / columns)
  body_top <- y + height - if (primary) 0.040 else 0.028
  body_bottom <- y + 0.010
  line_step <- min(0.0132, (body_top - body_bottom) / max(rows, 1L))
  col_width <- width / columns
  for (i in seq_along(functions)) {
    column <- (i - 1L) %/% rows
    row <- (i - 1L) %% rows
    item_x <- x + 0.008 + column * col_width
    item_y <- body_top - row * line_step
    if (primary) {
      grid::grid.circle(
        x = item_x, y = item_y, r = grid::unit(0.0018, "npc"),
        gp = grid::gpar(fill = colors$orange, col = colors$orange)
      )
    }
    grid::grid.text(
      paste0(functions[[i]], "()"),
      x = item_x + if (primary) 0.006 else 0,
      y = item_y, just = c("left", "center"),
      gp = grid::gpar(col = colors$ink, fontfamily = "mono", fontsize = 7.1)
    )
  }
  if (predict) {
    grid::grid.text(
      "predict()  [S3]", x = x + 0.014,
      y = body_top - rows * line_step, just = c("left", "center"),
      gp = grid::gpar(
        col = colors$orange, fontfamily = "mono", fontface = "bold", fontsize = 7.2
      )
    )
  }
}

draw_arrow <- function(x0, y0, x1, y1) {
  grid::grid.lines(
    x = c(x0, x1), y = c(y0, y1),
    arrow = grid::arrow(type = "closed", length = grid::unit(0.07, "inches")),
    gp = grid::gpar(col = colors$orange, lwd = 2.0, lineend = "round")
  )
}

draw_map <- function() {
  grid::grid.newpage()
  grid::grid.rect(gp = grid::gpar(fill = colors$background, col = NA))
  grid::grid.text(
    "rGEDI function architecture", x = 0.025, y = 0.975,
    just = c("left", "top"),
    gp = grid::gpar(col = colors$green, fontfamily = "sans", fontface = "bold", fontsize = 26)
  )
  grid::grid.text(
    sprintf("Complete public API: %d exported functions, grouped by workflow and support role", length(exports)),
    x = 0.025, y = 0.946, just = c("left", "top"),
    gp = grid::gpar(col = colors$muted, fontfamily = "sans", fontsize = 11)
  )

  xs <- c(0.025, 0.2675, 0.5100, 0.7525)
  w <- 0.2225
  top_y <- 0.765
  top_h <- 0.155
  for (i in 1:4) {
    draw_card(names(main_groups)[[i]], main_groups[[i]], xs[[i]], top_y, w, top_h)
  }
  for (i in 1:3) draw_arrow(xs[[i]] + w + 0.002, top_y + top_h / 2, xs[[i + 1L]] - 0.002, top_y + top_h / 2)

  lower_y <- 0.565
  lower_h <- 0.178
  draw_card(names(main_groups)[[5L]], main_groups[[5L]], 0.025, lower_y, 0.465, lower_h, columns = 2L)
  draw_card(names(main_groups)[[6L]], main_groups[[6L]], 0.510, lower_y, w, lower_h, predict = TRUE)
  draw_card(names(main_groups)[[7L]], main_groups[[7L]], 0.7525, lower_y, w, lower_h)
  draw_arrow(0.490, lower_y + lower_h / 2, 0.508, lower_y + lower_h / 2)
  draw_arrow(0.7325, lower_y + lower_h / 2, 0.7505, lower_y + lower_h / 2)
  grid::grid.lines(
    x = c(0.975, 0.988, 0.988, 0.012, 0.012, 0.023),
    y = c(top_y + top_h / 2, top_y + top_h / 2, 0.753, 0.753, lower_y + lower_h / 2, lower_y + lower_h / 2),
    arrow = grid::arrow(type = "closed", length = grid::unit(0.07, "inches")),
    gp = grid::gpar(col = colors$orange, lwd = 1.8, lineend = "round", linejoin = "round")
  )

  grid::grid.text(
    "Auxiliary functions", x = 0.025, y = 0.528,
    just = c("left", "center"),
    gp = grid::gpar(col = colors$green, fontfamily = "sans", fontface = "bold", fontsize = 15)
  )
  grid::grid.text(
    "Every exported helper is listed once", x = 0.170, y = 0.528,
    just = c("left", "center"),
    gp = grid::gpar(col = colors$muted, fontfamily = "sans", fontsize = 8)
  )

  aux_layout <- list(
    c("Sampling specifications", "Cloud and package configuration", "Classes and storage helpers"),
    c("Raster product operations", "Archive and batch processing"),
    c("GEE catalog, display and export", "GEE predictors and geometry")
  )
  aux_x <- c(0.025, 0.345, 0.665)
  aux_w <- 0.300
  for (column in seq_along(aux_layout)) {
    y_top <- 0.505
    for (group_name in aux_layout[[column]]) {
      values <- auxiliary_groups[[group_name]]
      card_h <- 0.040 + 0.0132 * length(values)
      y <- y_top - card_h
      draw_card(group_name, values, aux_x[[column]], y, aux_w, card_h, primary = FALSE)
      y_top <- y - 0.010
    }
  }

  s3_y <- 0.035
  s3_h <- 0.105
  draw_box(0.025, s3_y, 0.950, s3_h, colors$white, colors$line, lwd = 1.1)
  grid::grid.text(
    "Registered S3 methods", x = 0.034, y = s3_y + s3_h - 0.015,
    just = c("left", "center"),
    gp = grid::gpar(col = colors$green, fontfamily = "sans", fontface = "bold", fontsize = 8.5)
  )
  grid::grid.text(
    "Automatic methods used by collections, cloud references, models and printing",
    x = 0.145, y = s3_y + s3_h - 0.015, just = c("left", "center"),
    gp = grid::gpar(col = colors$muted, fontfamily = "sans", fontsize = 7.3)
  )
  method_rows <- 5L
  method_col_width <- 0.31
  for (i in seq_along(s3_methods)) {
    column <- (i - 1L) %/% method_rows
    row <- (i - 1L) %% method_rows
    grid::grid.text(
      s3_methods[[i]], x = 0.034 + column * method_col_width,
      y = s3_y + s3_h - 0.038 - row * 0.013,
      just = c("left", "center"),
      gp = grid::gpar(col = colors$ink, fontfamily = "mono", fontsize = 7.0)
    )
  }
  grid::grid.text(
    "Primary functions form the user workflow; source-specific helpers remain available for explicit control.",
    x = 0.025, y = 0.012, just = c("left", "center"),
    gp = grid::gpar(col = colors$orange, fontfamily = "sans", fontsize = 8)
  )
}

png_file <- file.path(output_dir, "rGEDI-function-map.png")
pdf_file <- file.path(output_dir, "rGEDI-function-map.pdf")

grDevices::png(png_file, width = 3200, height = 3580, res = 200, type = "cairo")
draw_map()
grDevices::dev.off()

grDevices::pdf(pdf_file, width = 16, height = 17.9, useDingbats = FALSE, onefile = TRUE)
draw_map()
grDevices::dev.off()

message(
  "Wrote ", png_file, " and ", pdf_file,
  "; verified ", length(exports), " exports exactly once."
)

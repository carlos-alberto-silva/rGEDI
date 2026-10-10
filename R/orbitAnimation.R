#' Animate a GEDI ground track around the Earth
#'
#' Creates an interactive 3D HTML animation or an animated GIF from GEDI
#' footprint coordinates. The presentation follows the ICESat2VegR orbit
#' animation, with a rotating Earth, a NASA ISS illustration carrying GEDI,
#' and a red near-infrared laser and ground track. ISS motion follows the
#' time-ordered geolocation of the most complete beam extracted from the
#' supplied HDF5 object.
#'
#' @param x A data frame, `data.table`, `sf` or `SpatVector` containing GEDI
#'   footprint coordinates, or an open `gedi.level1b`, `gedi.level2a`,
#'   `gedi.level2b`, or `gedi.level4a` object.
#' @param lon,lat Names of longitude and latitude columns. When `NULL`, common
#'   GEDI coordinate names are detected automatically.
#' @param time Optional name of a column used to order footprints within each
#'   track. GEDI `delta_time` and `shot_number` columns are detected by default.
#' @param track Optional name of a beam or track column. Separate tracks are
#'   colored independently.
#' @param output_file Path of the `.html` or `.gif` file to create. GIF output
#'   requires the suggested `gifski` package.
#' @param title Title shown above the animation.
#' @param duration Duration of one animation cycle in seconds.
#' @param launch Open the animation in the default browser after it is written.
#' @param track_speed Initial HTML track playback speed from 1 to 15.
#' @param earth_rotation_speed Initial Earth rotation speed from 0 to 20.
#'
#' @return The normalized path to the HTML file, invisibly.
#' @export
#'
#' @examples
#' orbit <- data.frame(
#'   longitude = seq(-60, -45, length.out = 30),
#'   latitude = seq(-8, 5, length.out = 30),
#'   delta_time = seq_len(30), beam = "BEAM0101"
#' )
#' html <- plot_gedi_orbit_animation(
#'   orbit, output_file = tempfile(fileext = ".html")
#' )
#' file.exists(html)
plot_gedi_orbit_animation <- function(
    x,
    lon = NULL,
    lat = NULL,
    time = NULL,
    track = NULL,
    output_file = tempfile("rGEDI-orbit-", fileext = ".html"),
    title = "GEDI orbit animation",
    duration = 18,
    launch = interactive(),
    track_speed = 2,
    earth_rotation_speed = 2) {
  duration <- as.numeric(duration)[1L]
  if (!is.finite(duration) || duration <= 0) {
    stop("`duration` must be a positive number of seconds.", call. = FALSE)
  }
  track_speed <- as.numeric(track_speed)[1L]
  earth_rotation_speed <- as.numeric(earth_rotation_speed)[1L]
  if (!is.finite(track_speed) || track_speed < 1 || track_speed > 15) {
    stop("`track_speed` must be between 1 and 15.", call. = FALSE)
  }
  if (!is.finite(earth_rotation_speed) || earth_rotation_speed < 0 ||
      earth_rotation_speed > 20) {
    stop("`earth_rotation_speed` must be between 0 and 20.", call. = FALSE)
  }

  output <- as.data.frame(getGEDITrack(
    x, lon = lon, lat = lat, time = time, track = track
  ))
  names(output)[match(c("longitude", "latitude", "sequence"), names(output))] <-
    c("lon", "lat", "order")
  reference_track <- names(which.max(table(output$track)))[1L]
  output$reference <- output$track == reference_track
  output$index <- seq_len(nrow(output))

  extension <- tolower(tools::file_ext(output_file))
  if (identical(extension, "gif")) {
    return(.write_gedi_orbit_gif(
      output, output_file, title, duration, launch, earth_rotation_speed
    ))
  }
  if (!identical(extension, "html")) {
    stop("`output_file` must end in .html or .gif.", call. = FALSE)
  }

  payload <- jsonlite::toJSON(
    output[c("lon", "lat", "track", "reference", "index")],
    dataframe = "rows", auto_unbox = TRUE, na = "null", digits = 8
  )
  payload <- gsub("</", "<\\\\/", payload, fixed = TRUE)
  safe_title <- .html_escape(title)
  html <- .orbit_html(
    payload, safe_title, duration, track_speed, earth_rotation_speed
  )
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(html, output_file, useBytes = TRUE)
  result <- normalizePath(output_file, winslash = "/", mustWork = TRUE)
  if (isTRUE(launch)) {
    utils::browseURL(result)
  }
  invisible(result)
}

#' Extract a GEDI reference ground track
#'
#' Standardizes footprint coordinates, acquisition order, and beam identifiers
#' from any open point-level GEDI product or an extracted footprint table. The
#' result can be plotted directly or passed to [plot_gedi_orbit_animation()].
#'
#' @inheritParams plot_gedi_orbit_animation
#' @param beams Optional beam names to read from an open GEDI object.
#' @param every Keep every nth footprint after ordering. This is useful when
#'   plotting a complete orbit.
#' @return A [data.table::data.table] with `longitude`, `latitude`, `sequence`,
#'   `track`, and available time and shot identifiers.
#' @export
#' @examples
#' shots <- data.frame(
#'   lon_lowestmode = seq(-44.2, -44.1, length.out = 10),
#'   lat_lowestmode = seq(-13.8, -13.7, length.out = 10),
#'   delta_time = seq_len(10), beam = "BEAM0000"
#' )
#' getGEDITrack(shots, every = 2)
getGEDITrack <- function(x, lon = NULL, lat = NULL, time = NULL, track = NULL,
                         beams = NULL, every = 1L) {
  points <- .gedi_orbit_points(x, beams = beams)
  choices <- names(points)
  lon <- lon %||orbit% .first_orbit_name(
    choices, c("longitude", "lon_lowestmode", "longitude_bin0", "lon")
  )
  lat <- lat %||orbit% .first_orbit_name(
    choices, c("latitude", "lat_lowestmode", "latitude_bin0", "lat")
  )
  if (is.null(lon) || is.null(lat) || !all(c(lon, lat) %in% choices)) {
    stop("Could not identify longitude and latitude columns.", call. = FALSE)
  }
  time <- time %||orbit% .first_orbit_name(
    choices, c("delta_time", "shot_number", "time", "datetime")
  )
  track <- track %||orbit% .first_orbit_name(
    choices, c("beam", "track", "track_id", "orbit")
  )
  every <- as.integer(every)[1L]
  if (!is.finite(every) || every < 1L) stop("`every` must be a positive integer.")
  keep <- is.finite(suppressWarnings(as.numeric(points[[lon]]))) &
    is.finite(suppressWarnings(as.numeric(points[[lat]])))
  points <- points[keep, , drop = FALSE]
  if (!nrow(points)) stop("No finite GEDI coordinates were supplied.", call. = FALSE)
  order_value <- if (is.null(time)) seq_len(nrow(points)) else points[[time]]
  track_value <- if (is.null(track)) rep("GEDI", nrow(points)) else as.character(points[[track]])
  track_value[is.na(track_value) | !nzchar(track_value)] <- "GEDI"
  index <- order(track_value, order_value, na.last = TRUE)
  index <- index[seq.int(1L, length(index), by = every)]
  ans <- data.table::data.table(
    longitude = as.numeric(points[[lon]][index]),
    latitude = as.numeric(points[[lat]][index]),
    sequence = order_value[index],
    track = track_value[index]
  )
  for (field in intersect(c("delta_time", "shot_number", "beam"), choices)) {
    if (!field %in% names(ans)) ans[[field]] <- points[[field]][index]
  }
  # A granule can contain separated passes with the same beam identifier. Split
  # large spatial jumps so plotting does not draw artificial cross-track lines.
  for (label in unique(ans$track)) {
    rows <- which(ans$track == label)
    if (length(rows) < 3L) next
    distance <- sqrt(diff(ans$longitude[rows])^2 + diff(ans$latitude[rows])^2)
    local_step <- stats::median(distance[is.finite(distance) & distance > 0], na.rm = TRUE)
    if (!is.finite(local_step) || local_step <= 0) next
    segment <- cumsum(c(TRUE, distance > 5 * local_step))
    if (max(segment) > 1L) ans$track[rows] <- paste0(label, ".", segment)
  }
  ans
}

`%||orbit%` <- function(x, y) if (is.null(x)) y else x

.first_orbit_name <- function(nms, candidates) {
  found <- candidates[candidates %in% nms]
  if (length(found)) found[[1L]] else NULL
}

.gedi_orbit_points <- function(x, beams = NULL) {
  if (methods::is(x, "gedi.level1b")) {
    return(as.data.frame(getLevel1BGeo(x, select = "delta_time", beams = beams)))
  }
  if (methods::is(x, "gedi.level2a")) {
    return(as.data.frame(getLevel2AM(x, beams = beams, include_rh = FALSE)))
  }
  if (methods::is(x, "gedi.level2b")) {
    return(as.data.frame(getLevel2BVPM(
      x, cols = c("beam", "shot_number", "delta_time",
                  "latitude_bin0", "longitude_bin0"), beams = beams
    )))
  }
  if (methods::is(x, "gedi.level4a")) {
    return(as.data.frame(getLevel4A(
      x, cols = c("shot_number", "delta_time", "lat_lowestmode",
                  "lon_lowestmode"), beams = beams
    )))
  }
  if (inherits(x, "SpatVector")) {
    crds <- terra::crds(terra::project(x, "EPSG:4326"))
    return(cbind(as.data.frame(x), longitude = crds[, 1L], latitude = crds[, 2L]))
  }
  if (inherits(x, "sf")) {
    transformed <- sf::st_transform(x, 4326)
    crds <- sf::st_coordinates(transformed)
    output <- sf::st_drop_geometry(transformed)
    output$longitude <- crds[, 1L]
    output$latitude <- crds[, 2L]
    return(as.data.frame(output))
  }
  if (!is.data.frame(x)) {
    stop("`x` must contain GEDI coordinates or be an open GEDI object.", call. = FALSE)
  }
  as.data.frame(x)
}

.orbit_asset <- function(name) {
  installed <- system.file("extdata", "orbit", name, package = "rGEDI")
  if (nzchar(installed) && file.exists(installed)) return(installed)
  source <- file.path("inst", "extdata", "orbit", name)
  if (file.exists(source)) return(normalizePath(source, mustWork = TRUE))
  stop("Orbit animation asset not found: ", name, call. = FALSE)
}

.orbit_data_uri <- function(path, mime) {
  bytes <- readBin(path, what = "raw", n = file.info(path)$size)
  encoded <- jsonlite::base64_enc(bytes)
  encoded <- paste(encoded, collapse = "")
  encoded <- gsub(intToUtf8(10), "", encoded, fixed = TRUE)
  encoded <- gsub(intToUtf8(13), "", encoded, fixed = TRUE)
  paste0("data:", mime, ";base64,", encoded)
}

.write_gedi_orbit_gif <- function(output, output_file, title, duration, launch,
                                  earth_rotation_speed = 2) {
  if (!requireNamespace("gifski", quietly = TRUE)) {
    stop("GIF output requires the suggested 'gifski' package.", call. = FALSE)
  }
  orbit <- output[output$reference, , drop = FALSE]
  orbit <- orbit[order(orbit$order, na.last = TRUE), , drop = FALSE]
  nframes <- min(48L, nrow(orbit))
  frame_rows <- unique(round(seq(1, nrow(orbit), length.out = nframes)))
  frames <- file.path(tempdir(), sprintf("rGEDI-orbit-%03d.png", seq_along(frame_rows)))
  on.exit(unlink(frames, force = TRUE), add = TRUE)
  colors <- stats::setNames(rep("#ff1744", length(unique(output$track))),
                            unique(output$track))
  center_lon <- stats::median(output$lon, na.rm = TRUE)
  center_lat <- max(-35, min(35, stats::median(output$lat, na.rm = TRUE)))
  world <- if (requireNamespace("maps", quietly = TRUE)) {
    maps::map("world", plot = FALSE, fill = TRUE)
  } else NULL
  set.seed(42)
  stars <- data.frame(x = runif(180, -1.55, 1.55), y = runif(180, -1.1, 1.1))
  iss_image <- if (requireNamespace("png", quietly = TRUE)) {
    png::readPNG(.orbit_asset("ISS-NASA-transparent.png"))
  } else NULL
  for (i in seq_along(frame_rows)) {
    row <- frame_rows[[i]]
    progress <- if (nrow(orbit) == 1L) 1 else (row - 1) / (nrow(orbit) - 1)
    rotation <- center_lon + 360 * progress * earth_rotation_speed / 2
    grDevices::png(frames[[i]], width = 1000, height = 650, bg = "#02060b")
    graphics::par(mar = c(0, 0, 2.2, 0), fg = "white")
    graphics::plot.new()
    graphics::plot.window(c(-1.58, 1.58), c(-1.08, 1.08), asp = 1)
    graphics::points(stars$x, stars$y, pch = ".", col = "#ffffff99")
    graphics::symbols(0, 0, circles = 1, inches = FALSE, add = TRUE,
                      bg = "#2387a4", fg = "#9edce8", lwd = 2)
    .draw_orbit_graticule(rotation, center_lat)
    if (!is.null(world)) .draw_orbit_land(world, rotation, center_lat)
    for (beam in unique(output$track)) {
      part <- output[output$track == beam, , drop = FALSE]
      part <- part[order(part$order, na.last = TRUE), , drop = FALSE]
      part <- part[seq_len(max(1L, round(progress * nrow(part)))), , drop = FALSE]
      projected <- .orbit_project(part$lon, part$lat, rotation, center_lat, 1.018)
      .draw_visible_orbit_line(projected, colors[[beam]], 3)
    }
    focus_row <- orbit[row, , drop = FALSE]
    focus <- .orbit_project(focus_row$lon, focus_row$lat,
                            rotation, center_lat, 1.12)
    ground <- .orbit_project(focus_row$lon, focus_row$lat,
                             rotation, center_lat, 1)
    if (isTRUE(focus$visible)) {
      graphics::segments(focus$x, focus$y, ground$x, ground$y,
                         col = "#ff1744", lwd = 2)
      if (is.null(iss_image)) {
        .draw_iss(focus$x, focus$y)
      } else {
        graphics::rasterImage(
          iss_image, focus$x - 0.15, focus$y - 0.09,
          focus$x + 0.15, focus$y + 0.09, interpolate = TRUE
        )
      }
      graphics::text(focus$x, focus$y + 0.12, "ISS + GEDI", col = "white",
                     cex = 0.9, font = 2)
    }
    graphics::title(main = title, col.main = "white", cex.main = 1.35)
    graphics::text(1.48, -1.01,
                   sprintf("HDF5 orbit point %s of %s", row, nrow(orbit)),
                   adj = 1, col = "#b6d6df", cex = 0.82)
    grDevices::dev.off()
  }
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  gifski::gifski(
    frames, gif_file = output_file, width = 900, height = 600,
    delay = duration / length(frames), loop = TRUE, progress = FALSE
  )
  result <- normalizePath(output_file, winslash = "/", mustWork = TRUE)
  if (isTRUE(launch)) utils::browseURL(result)
  invisible(result)
}

.orbit_project <- function(lon, lat, center_lon, center_lat = 0, radius = 1) {
  lambda <- (as.numeric(lon) - center_lon) * pi / 180
  phi <- as.numeric(lat) * pi / 180
  phi0 <- center_lat * pi / 180
  visible <- sin(phi0) * sin(phi) + cos(phi0) * cos(phi) * cos(lambda) >= 0
  list(
    x = radius * cos(phi) * sin(lambda),
    y = radius * (cos(phi0) * sin(phi) - sin(phi0) * cos(phi) * cos(lambda)),
    visible = visible
  )
}

.draw_visible_orbit_line <- function(projected, col, lwd = 1) {
  group <- cumsum(c(TRUE, diff(as.integer(projected$visible)) != 0))
  for (g in unique(group[projected$visible])) {
    keep <- group == g & projected$visible
    if (sum(keep) > 1) graphics::lines(projected$x[keep], projected$y[keep],
                                       col = col, lwd = lwd)
  }
}

.draw_orbit_graticule <- function(center_lon, center_lat) {
  for (latitude in seq(-60, 60, 30)) {
    p <- .orbit_project(seq(-180, 180, 2), latitude, center_lon, center_lat)
    .draw_visible_orbit_line(p, "#d7f5fa33", 0.8)
  }
  for (longitude in seq(-180, 150, 30)) {
    p <- .orbit_project(longitude, seq(-90, 90, 2), center_lon, center_lat)
    .draw_visible_orbit_line(p, "#d7f5fa33", 0.8)
  }
}

.draw_orbit_land <- function(world, center_lon, center_lat) {
  breaks <- c(0L, which(is.na(world$x)), length(world$x) + 1L)
  for (i in seq_len(length(breaks) - 1L)) {
    rows <- seq.int(breaks[i] + 1L, breaks[i + 1L] - 1L)
    if (length(rows) < 3L) next
    p <- .orbit_project(world$x[rows], world$y[rows], center_lon, center_lat)
    visible <- which(p$visible)
    if (length(visible) >= 3L) {
      graphics::polygon(p$x[visible], p$y[visible], col = "#5e8f56",
                        border = "#b7d49b88", lwd = 0.4)
    }
  }
}

.draw_iss <- function(x, y) {
  # A compact ISS silhouette: two solar arrays, central truss, and GEDI payload.
  graphics::rect(x - 0.09, y - 0.018, x - 0.025, y + 0.018,
                 col = "#2878b8", border = "white", lwd = 0.7)
  graphics::rect(x + 0.025, y - 0.018, x + 0.09, y + 0.018,
                 col = "#2878b8", border = "white", lwd = 0.7)
  graphics::segments(x - 0.105, y, x + 0.105, y, col = "#e8edf0", lwd = 2)
  graphics::rect(x - 0.023, y - 0.024, x + 0.023, y + 0.024,
                 col = "#e8edf0", border = "#25333d")
  graphics::points(x, y - 0.028, pch = 21, bg = "#f2b134", col = "white", cex = 0.8)
}

.html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  gsub("'", "&#39;", x, fixed = TRUE)
}

.orbit_html <- function(payload, title, duration, track_speed = 2,
                        earth_rotation_speed = 2) {
  template <- .orbit_template_three()
  template <- sub("__GEDI_DATA__", payload, template, fixed = TRUE)
  template <- sub("__GEDI_TITLE__", title, template, fixed = TRUE)
  template <- sub("__GEDI_DURATION__", format(duration, scientific = FALSE),
                  template, fixed = TRUE)
  template <- sub("__TRACK_SPEED__", format(track_speed, scientific = FALSE),
                  template, fixed = TRUE)
  template <- sub("__ROTATION_SPEED__",
                  format(earth_rotation_speed, scientific = FALSE), template,
                  fixed = TRUE)
  earth_uri <- .orbit_data_uri(
    .orbit_asset("Stylized_World_Topo_5400x2700.jpeg"), "image/jpeg"
  )
  iss_uri <- .orbit_data_uri(
    .orbit_asset("ISS-NASA-transparent.png"), "image/png"
  )
  template <- sub("__EARTH_TEXTURE__", earth_uri, template, fixed = TRUE)
  sub("__ISS_IMAGE__", iss_uri, template, fixed = TRUE)
}

.orbit_template_three <- function() '<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__GEDI_TITLE__</title>
<style>
html,body{margin:0;padding:0;width:100%;height:100%;overflow:hidden;background:#000;font-family:Arial,sans-serif}
canvas{display:block}#controls{position:absolute;top:15px;left:15px;z-index:10;background:rgba(0,0,0,.78);color:#fff;padding:14px;border-radius:10px;width:450px;max-width:calc(100vw - 58px);font-size:14px;border:1px solid #ffffff30}
button{margin:3px;padding:6px 10px;border:0;border-radius:5px;cursor:pointer}input[type=range]{width:330px;max-width:70vw;accent-color:#ff1744}
#timeLabel,#trackLabel,#status{margin-top:8px}.credit{margin-top:10px;color:#b7c9d8;font-size:11px}
</style></head><body>
<div id="controls"><b>__GEDI_TITLE__</b><br><br>
<button onclick="playAnimation()">Play</button><button onclick="pauseAnimation()">Pause</button><button onclick="resetAnimation()">Reset</button><br><br>
Track speed:<br><input type="range" min="1" max="15" value="__TRACK_SPEED__" id="trackSpeed"> <span id="trackSpeedValue">__TRACK_SPEED__</span><br><br>
Earth rotation speed:<br><input type="range" min="0" max="20" value="__ROTATION_SPEED__" id="rotationSpeed"> <span id="rotationSpeedValue">__ROTATION_SPEED__</span>
<div id="trackLabel">Track:</div><div id="timeLabel">Coordinates:</div><div id="status">Loading...</div>
<div class="credit">ISS illustration: NASA | GEDI wavelength: 1064 nm near infrared</div></div>
<script src="https://cdn.jsdelivr.net/npm/three@0.128.0/build/three.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
<script>
const allGEDIData=__GEDI_DATA__;
const trackData=allGEDIData.filter(function(d){return d.reference;}).sort(function(a,b){return a.index-b.index;});
const earthTexture="__EARTH_TEXTURE__",issImage="__ISS_IMAGE__";
const scene=new THREE.Scene();scene.background=new THREE.Color(0x000000);
const camera=new THREE.PerspectiveCamera(45,window.innerWidth/window.innerHeight,.1,1000);camera.position.set(0,0,5);
const renderer=new THREE.WebGLRenderer({antialias:true,alpha:false});renderer.setSize(window.innerWidth,window.innerHeight);renderer.setPixelRatio(window.devicePixelRatio);document.body.appendChild(renderer.domElement);
const controls=new THREE.OrbitControls(camera,renderer.domElement);controls.enableDamping=true;controls.dampingFactor=.05;controls.enableZoom=true;
scene.add(new THREE.AmbientLight(0xffffff,2));
const earthGroup=new THREE.Group();scene.add(earthGroup);const earthRadius=1.5,orbitRadius=earthRadius+.35;
const earth=new THREE.Mesh(new THREE.SphereGeometry(earthRadius,256,256),new THREE.MeshBasicMaterial({color:0x1e66b1}));earthGroup.add(earth);
new THREE.TextureLoader().load(earthTexture,function(texture){texture.minFilter=THREE.LinearFilter;texture.magFilter=THREE.LinearFilter;texture.generateMipmaps=false;earth.material.dispose();earth.material=new THREE.MeshBasicMaterial({map:texture});document.getElementById("status").innerHTML="Earth texture and GEDI orbit loaded";});
const starGeometry=new THREE.BufferGeometry(),starPositions=[];for(let i=0;i<7000;i++){starPositions.push((Math.random()-.5)*120,(Math.random()-.5)*120,(Math.random()-.5)*120);}starGeometry.setAttribute("position",new THREE.Float32BufferAttribute(starPositions,3));scene.add(new THREE.Points(starGeometry,new THREE.PointsMaterial({color:0xffffff,size:.04})));
function latLonToVector3(lat,lon,radius){const phi=(90-lat)*Math.PI/180,theta=(lon+180)*Math.PI/180;return new THREE.Vector3(-radius*Math.sin(phi)*Math.cos(theta),radius*Math.cos(phi),radius*Math.sin(phi)*Math.sin(theta));}
const issTexture=new THREE.TextureLoader().load(issImage),satellite=new THREE.Sprite(new THREE.SpriteMaterial({map:issTexture,transparent:true,depthWrite:false}));satellite.scale.set(.58,.34,1);earthGroup.add(satellite);
const labelCanvas=document.createElement("canvas");labelCanvas.width=1024;labelCanvas.height=256;const labelContext=labelCanvas.getContext("2d");labelContext.font="bold 90px Arial";labelContext.textAlign="center";labelContext.textBaseline="middle";labelContext.lineWidth=10;labelContext.strokeStyle="black";labelContext.strokeText("ISS + GEDI",512,128);labelContext.fillStyle="white";labelContext.fillText("ISS + GEDI",512,128);
const satelliteLabel=new THREE.Sprite(new THREE.SpriteMaterial({map:new THREE.CanvasTexture(labelCanvas),transparent:true,depthWrite:false,depthTest:false}));satelliteLabel.scale.set(.8,.2,1);earthGroup.add(satelliteLabel);
let laserBeam=null,activeTrackLine=null,pointIndex=0,running=true,trackPoints=[];const trackColor=0xff1744;
function updateLaserBeam(satPos){if(laserBeam){earthGroup.remove(laserBeam);laserBeam.geometry.dispose();laserBeam.material.dispose();}const ground=satPos.clone().normalize().multiplyScalar(earthRadius+.01);laserBeam=new THREE.Line(new THREE.BufferGeometry().setFromPoints([satPos,ground]),new THREE.LineBasicMaterial({color:0xff1744,transparent:true,opacity:1}));earthGroup.add(laserBeam);}
function clearTrackLine(){if(activeTrackLine){earthGroup.remove(activeTrackLine);activeTrackLine.geometry.dispose();activeTrackLine.material.dispose();activeTrackLine=null;}}
function updateTrackLine(){clearTrackLine();if(trackPoints.length<2)return;activeTrackLine=new THREE.Line(new THREE.BufferGeometry().setFromPoints(trackPoints),new THREE.LineBasicMaterial({color:trackColor}));earthGroup.add(activeTrackLine);}
function updateSatellitePosition(pos){satellite.position.copy(pos);satelliteLabel.position.copy(pos.clone().add(pos.clone().normalize().multiplyScalar(.28)));updateLaserBeam(pos);}
function addCurrentPoint(){if(pointIndex>=trackData.length){running=false;document.getElementById("status").innerHTML="GEDI track completed";return;}const p=trackData[pointIndex],pos=latLonToVector3(p.lat,p.lon,orbitRadius);trackPoints.push(pos);updateSatellitePosition(pos);document.getElementById("trackLabel").innerHTML="HDF5 orbit point: "+(pointIndex+1)+" of "+trackData.length+" | Beam: "+p.track;document.getElementById("timeLabel").innerHTML="Longitude: "+p.lon.toFixed(4)+"&deg; | Latitude: "+p.lat.toFixed(4)+"&deg;";pointIndex++;updateTrackLine();}
function playAnimation(){running=true;document.getElementById("status").innerHTML="Animation running";}
function pauseAnimation(){running=false;document.getElementById("status").innerHTML="Animation paused";}
function resetAnimation(){running=false;pointIndex=0;trackPoints=[];clearTrackLine();if(laserBeam){earthGroup.remove(laserBeam);laserBeam.geometry.dispose();laserBeam.material.dispose();laserBeam=null;}if(trackData.length)addCurrentPoint();document.getElementById("status").innerHTML="Animation reset";}
document.getElementById("trackSpeed").addEventListener("input",function(){document.getElementById("trackSpeedValue").innerHTML=this.value;});document.getElementById("rotationSpeed").addEventListener("input",function(){document.getElementById("rotationSpeedValue").innerHTML=this.value;});if(trackData.length)addCurrentPoint();else document.getElementById("status").innerHTML="No GEDI orbit points found";
function animate(){requestAnimationFrame(animate);const rotationSpeed=Number(document.getElementById("rotationSpeed").value);earthGroup.rotation.y+=.0008*rotationSpeed;if(running&&pointIndex<trackData.length){const trackSpeed=Number(document.getElementById("trackSpeed").value);for(let s=0;s<trackSpeed;s++){if(running&&pointIndex<trackData.length)addCurrentPoint();}}controls.update();renderer.render(scene,camera);}
animate();window.addEventListener("resize",function(){camera.aspect=window.innerWidth/window.innerHeight;camera.updateProjectionMatrix();renderer.setSize(window.innerWidth,window.innerHeight);});
</script></body></html>'

#' Animate a GEDI ground track around the Earth
#'
#' Creates a self-contained HTML animation or an animated GIF from GEDI
#' footprint coordinates. HTML output uses an interactive orthographic globe;
#' GIF output shows progressive reference-ground-track playback.
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
    launch = interactive()) {
  duration <- as.numeric(duration)[1L]
  if (!is.finite(duration) || duration <= 0) {
    stop("`duration` must be a positive number of seconds.", call. = FALSE)
  }

  output <- as.data.frame(getGEDITrack(
    x, lon = lon, lat = lat, time = time, track = track
  ))
  names(output)[match(c("longitude", "latitude", "sequence"), names(output))] <-
    c("lon", "lat", "order")
  output$index <- seq_len(nrow(output))

  extension <- tolower(tools::file_ext(output_file))
  if (identical(extension, "gif")) {
    return(.write_gedi_orbit_gif(output, output_file, title, duration, launch))
  }
  if (!identical(extension, "html")) {
    stop("`output_file` must end in .html or .gif.", call. = FALSE)
  }

  payload <- jsonlite::toJSON(
    output[c("lon", "lat", "track", "index")],
    dataframe = "rows", auto_unbox = TRUE, na = "null", digits = 8
  )
  payload <- gsub("</", "<\\\\/", payload, fixed = TRUE)
  safe_title <- .html_escape(title)
  html <- .orbit_html(payload, safe_title, duration)
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
    return(as.data.frame(getLevel2AM(x, beams = beams)))
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

.write_gedi_orbit_gif <- function(output, output_file, title, duration, launch) {
  if (!requireNamespace("gifski", quietly = TRUE)) {
    stop("GIF output requires the suggested 'gifski' package.", call. = FALSE)
  }
  nframes <- min(80L, nrow(output))
  frame_rows <- unique(round(seq(1, nrow(output), length.out = nframes)))
  frames <- file.path(tempdir(), sprintf("rGEDI-orbit-%03d.png", seq_along(frame_rows)))
  on.exit(unlink(frames, force = TRUE), add = TRUE)
  xlim <- range(output$lon, finite = TRUE)
  ylim <- range(output$lat, finite = TRUE)
  colors <- grDevices::hcl.colors(length(unique(output$track)), "Dark 3")
  names(colors) <- unique(output$track)
  for (i in seq_along(frame_rows)) {
    row <- frame_rows[[i]]
    grDevices::png(frames[[i]], width = 900, height = 600, bg = "#07131d")
    graphics::par(mar = c(4.5, 4.8, 3.5, 1.5), fg = "white",
                  col.axis = "white", col.lab = "white", col.main = "white")
    graphics::plot(NA, xlim = xlim, ylim = ylim, asp = 1,
                   xlab = "Longitude", ylab = "Latitude", main = title)
    graphics::grid(col = "#ffffff28")
    shown <- output[seq_len(row), , drop = FALSE]
    for (beam in unique(shown$track)) {
      part <- shown[shown$track == beam, , drop = FALSE]
      graphics::lines(part$lon, part$lat, col = colors[[beam]], lwd = 2)
    }
    graphics::points(shown$lon[nrow(shown)], shown$lat[nrow(shown)],
                     pch = 21, bg = "#ff5b45", col = "white", cex = 1.8)
    graphics::mtext(sprintf("Footprint %s of %s", row, nrow(output)),
                    side = 3, adj = 1, col = "#a9c5bd")
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

.html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  gsub("'", "&#39;", x, fixed = TRUE)
}

.orbit_html <- function(payload, title, duration) {
  template <- .orbit_template()
  template <- sub("__GEDI_DATA__", payload, template, fixed = TRUE)
  template <- sub("__GEDI_TITLE__", title, template, fixed = TRUE)
  sub("__GEDI_DURATION__", format(duration, scientific = FALSE), template,
      fixed = TRUE)
}

.orbit_template <- function() '<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__GEDI_TITLE__</title>
<style>
html,body{margin:0;height:100%;background:#07131d;color:#eef7f2;font-family:system-ui,sans-serif}
main{height:100%;display:grid;grid-template-rows:auto 1fr auto;overflow:hidden}
h1{font-size:clamp(18px,2.4vw,30px);font-weight:600;margin:18px 24px 4px}
p{margin:0 24px 12px;color:#a9c5bd}.stage{min-height:0;position:relative}
canvas{width:100%;height:100%;display:block}.panel{position:absolute;right:20px;top:10px;
background:#0b2230dd;border:1px solid #315568;border-radius:10px;padding:10px 14px;font-size:13px}
footer{display:flex;gap:12px;align-items:center;padding:10px 24px 18px}
button{background:#22b573;color:#fff;border:0;border-radius:7px;padding:7px 14px;cursor:pointer}
input{width:min(520px,60vw);accent-color:#22b573}
</style></head><body><main><header><h1>__GEDI_TITLE__</h1>
<p>Interactive GEDI ground-track playback</p></header><section class="stage">
<canvas id="globe"></canvas><div class="panel" id="readout"></div></section>
<footer><button id="play">Pause</button><input id="progress" type="range" min="0" max="1000" value="0"></footer>
</main><script>
const points=__GEDI_DATA__,duration=__GEDI_DURATION__*1000;
const canvas=document.getElementById("globe"),ctx=canvas.getContext("2d");
const slider=document.getElementById("progress"),button=document.getElementById("play"),readout=document.getElementById("readout");
let playing=true,start=performance.now(),manual=0;
const colors=["#51e2a5","#ffd166","#ef476f","#70d6ff","#c77dff","#ff9f1c","#90be6d","#f28482"];
const tracks=[...new Set(points.map(d=>d.track))],color=t=>colors[Math.max(0,tracks.indexOf(t))%colors.length];
function resize(){const d=devicePixelRatio||1,r=canvas.getBoundingClientRect();canvas.width=r.width*d;canvas.height=r.height*d;ctx.setTransform(d,0,0,d,0,0)}
addEventListener("resize",resize);resize();
function project(lon,lat,center,R,cx,cy){const a=(lon-center)*Math.PI/180,b=lat*Math.PI/180;
 const visible=Math.cos(b)*Math.cos(a)>=0;return{x:cx+R*Math.cos(b)*Math.sin(a),y:cy-R*Math.sin(b),visible};}
function draw(ts){const w=canvas.clientWidth,h=canvas.clientHeight,cx=w/2,cy=h/2,R=Math.max(40,Math.min(w,h)*.4);
 ctx.clearRect(0,0,w,h);const elapsed=playing?(ts-start)%duration:manual*duration;const p=elapsed/duration;
 if(playing)slider.value=Math.round(p*1000);const k=Math.min(points.length-1,Math.floor(p*points.length));
 const focus=points[k]||points[0],center=(focus?focus.lon:0)-25;
 const grad=ctx.createRadialGradient(cx-R*.3,cy-R*.4,R*.05,cx,cy,R);grad.addColorStop(0,"#2d8291");grad.addColorStop(.55,"#174e62");grad.addColorStop(1,"#071a29");
 ctx.beginPath();ctx.arc(cx,cy,R,0,Math.PI*2);ctx.fillStyle=grad;ctx.fill();ctx.save();ctx.beginPath();ctx.arc(cx,cy,R,0,Math.PI*2);ctx.clip();
 ctx.strokeStyle="#8bc5c52e";ctx.lineWidth=1;
 for(let lat=-60;lat<=60;lat+=30){ctx.beginPath();let pen=false;for(let lon=-180;lon<=180;lon+=3){const q=project(lon,lat,center,R,cx,cy);if(q.visible){pen?ctx.lineTo(q.x,q.y):ctx.moveTo(q.x,q.y);pen=true}else pen=false}ctx.stroke()}
 for(let lon=-180;lon<180;lon+=30){ctx.beginPath();let pen=false;for(let lat=-90;lat<=90;lat+=3){const q=project(lon,lat,center,R,cx,cy);if(q.visible){pen?ctx.lineTo(q.x,q.y):ctx.moveTo(q.x,q.y);pen=true}else pen=false}ctx.stroke()}
 for(const track of tracks){ctx.strokeStyle=color(track);ctx.lineWidth=2;ctx.beginPath();let pen=false;
  points.filter(d=>d.track===track).forEach(d=>{const q=project(d.lon,d.lat,center,R,cx,cy);if(q.visible){pen?ctx.lineTo(q.x,q.y):ctx.moveTo(q.x,q.y);pen=true}else pen=false});ctx.stroke()}
 ctx.restore();ctx.beginPath();ctx.arc(cx,cy,R,0,Math.PI*2);ctx.strokeStyle="#b4ece0aa";ctx.lineWidth=1.5;ctx.stroke();
 if(focus){const q=project(focus.lon,focus.lat,center,R,cx,cy);if(q.visible){ctx.beginPath();ctx.moveTo(q.x,q.y);ctx.lineTo(q.x,q.y-35);ctx.strokeStyle="#ff5b45aa";ctx.lineWidth=2;ctx.stroke();ctx.beginPath();ctx.arc(q.x,q.y-40,6,0,Math.PI*2);ctx.fillStyle="#fff";ctx.fill()}
 readout.innerHTML="Track: <b>"+focus.track+"</b><br>Longitude: "+focus.lon.toFixed(4)+"&deg;<br>Latitude: "+focus.lat.toFixed(4)+"&deg;"}
 requestAnimationFrame(draw)}
button.onclick=()=>{playing=!playing;button.textContent=playing?"Pause":"Play";if(playing)start=performance.now()-manual*duration};
slider.oninput=()=>{manual=slider.value/1000;playing=false;button.textContent="Play"};requestAnimationFrame(draw);
</script></body></html>'

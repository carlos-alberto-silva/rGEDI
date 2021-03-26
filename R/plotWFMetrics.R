#'GEDI full waveform plot with metrics
#'
#'@description Plots the waveform with overlaid RH metrics
#'
#'@usage plotWFMetrics(level1b, level2a, shot_number, rh=c(25, 50, 75),...)
#'
#'@param level1b A GEDI Level1B object (output of [readLevel1B()] function).
#'An S4 object of class "gedi.level1b".
#'@param level2a A GEDI Level2A object (output of [readLevel2A()] function).
#'An S4 object of class "gedi.level2a".
#'@param shot_number Shot number. A scalar representing the shot number of a giving pulse.
#'@param rh Integer vector. Specify which RH metrics to plot except rh0 and rh100, default c(25, 50, 75).
#'@param wf_args List. A list of arguments to pass to the waveform `graphics::plot()`.
#'@param cumEn_args List. A list of arguments to pass to the cumulative energy curve `graphics::plot()`.
#'@param height_axis_args List. A list of arguments to pass to the height `graphics::axis()`.
#'@param height_axis_text_args List. A list of arguments to pass to the height_axis label written as `graphics::text()`.
#'@param rh_arrow_args List. A list of arguments to pass to the rh `graphics::arrows()`.
#'@param rh_text_args List. A list of arguments to pass to the rh `graphics::text()`.
#'@param rh0_text_args List. A list of arguments to pass to the rh0 `graphics::text()`.
#'@param rh0_abline_args List. A list of arguments to pass to the rh0 `graphics::abline()`.
#'@param rh100_text_args List. A list of arguments to pass to the rh100 `graphics::text()`.
#'@param rh100_abline_args List. A list of arguments to pass to the rh100 `graphics::abline()`.
#'@param gz_text_args List. A list of arguments to pass to the gz `graphics::text()`.
#'@param gz_abline_args List. A list of arguments to pass to the gz `graphics::abline()`.
#'@param par_args List. A list of arguments to pass to the `graphics::par()`.
#'@param ... Will be passed to the main wf plot.
#'
#'@return Returns a raster layer(s) of selected GEDI Canopy Cover and Vertical Profile Metric(s)
#'
#'@seealso \url{https://lpdaac.usgs.gov/products/gedi02_bv002/}
#'
#'@examples
#'# specify the path to GEDI level1B and Level2A data (zip file)
#'outdir = tempdir()
#'level1B_fp_zip <- system.file("extdata",
#'                   "GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.zip",
#'                   package="rGEDI")
#'
#'level2A_fp_zip <- system.file("extdata",
#'                   "GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
#'                   package="rGEDI")
#'
#'# Unzipping GEDI level1B data
#'level1Bpath <- unzip(level1B_fp_zip,exdir = outdir)
#'level2Apath <- unzip(level2A_fp_zip,exdir = outdir)
#'
#'# Reading GEDI level1B and Level2A data (h5 file)
#'level1b<-readLevel1B(level1Bpath=level1Bpath)
#'level2a<-readLevel2A(level2Apath=level2Apath)
#'
#'shot_number = "19640521100108408"
#'
#'#
#'plotWFMetrics(level1b, level2a, shot_number,
#'              rh = c(25, 50, 75, 90),
#'              main = sprintf("Shot number: %s",
#'                           shot_number))
#'
#'# Only elevation axis
#'plotWFMetrics(level1b, level2a, shot_number,
#'              rh=c(25, 50, 75, 90),
#'              par_args = list(
#'                mar = c(5, 4, 4, 1)),
#'              height_axis_args = list(col=NA, col.ticks=NA),
#'              main=sprintf("Shot Number %s",
#'                           shot_number),
#')
#'#'
#'# Only height axis
#'plotWFMetrics(level1b, level2a, shot_number,
#'             rh=c(25, 50, 75, 90),
#'             wf_args = list(
#'               ylab=NA, yaxt="n",
#'               main=sprintf("Shot Number %s",
#'                           shot_number)),
#'             par_args = list(
#'               mar = c(5, 0, 4, 1),
#'               oma = c(0,4,0,0)),
#'             height_axis_args = list(side=2),
#'             height_axis_text_args = list(
#'               x=grconvertX(0, "npc"),
#'               adj=c(0.5, -4.7))
#')
#'
#'close(level1b)
#'close(level2a)
#'@importFrom stats quantile
#'@importFrom graphics abline arrows axis mtext par text
#'@export
plotWFMetrics = function(level1b, level2a, shot_number, rh=c(25, 50, 75),
                         wf_args=list(),
                         cumEn_args=list(),
                         height_axis_args=list(),
                         height_axis_text_args=list(),
                         rh_arrow_args=list(),
                         rh_text_args=list(),
                         rh0_text_args=list(),
                         rh0_abline_args=list(),
                         rh100_text_args=list(),
                         rh100_abline_args=list(),
                         gz_text_args=list(),
                         gz_abline_args=list(),
                         par_args=list(
                           mar = c(5, 4, 4, 5) + 0.3, oma=c(0, 0, 0, 0)
                         ),
                         height_nticks = 6,
                         ...) {
  # Avoid NOTEs from checking
  elevation = NULL
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  # Extracting GEDI full waveform for a giving shotnumber
  wf <- getLevel1BWF(level1b, shot_number=shot_number)

  level2AM<-getLevel2AM(level2a)
  shotid_mask = which(level2AM$shot_number == shot_number)

  level2_shot = level2AM[shotid_mask,]
  ground_z = level2_shot$elev_lowestmode
  rhs = paste0(as.character(c("rh0", as.character(paste0("rh",rh)), "rh100")))
  rh = level2_shot[, rhs, with=FALSE]
  rh_z = rh + ground_z


  top_z = level2AM[shotid_mask,]$elev_highestreturn

  range_energy = range(wf@dt$rxwaveform)
  range_abs_diff = abs(diff(range_energy))
  requireNamespace("data.table")

  range_z = range(rh_z)
  min_z = range_z[1]
  max_z = range_z[2]
  diff_z = max_z - min_z
  wf_between = wf@dt[elevation %between% range_z,,]
  energy_offset = min(range_energy)
  energy_no_offset = (wf_between$rxwaveform - energy_offset)
  cumsum_energy = cumsum(rev(energy_no_offset))

  range_cumsum = range(cumsum_energy)
  range_abs_diff_cumsum = abs(diff(range_cumsum))
  energy_cum_normalized = ((cumsum_energy)/(range_abs_diff_cumsum/range_abs_diff))+energy_offset

  par_default = list(mar = c(5, 4, 4, 5) + 0.3, oma=c(0, 0, 0, 0))
  par_default = modifyList(par_default, par_args)
  do.call(par, par_default)

  offset = diff_z*0.2
  ymin = min_z-offset
  ymax = max_z+offset
  wf_interest=wf@dt[wf@dt$elevation >= ymin & wf@dt$elevation <= ymax,]$rxwaveform
  qts=quantile(wf_interest, c(0.05, 1), type=1)

  z_masked = rev(wf_between$elevation)

  ticks = seq(min_z, max_z, length=height_nticks)
  closest_to_ground = which.min(abs(ticks-ground_z))
  ticks[closest_to_ground] = ground_z
  ticks_label = format(ticks-ground_z, digits = 2)


  rh_closest_en = list()
  for (i in 1:length(rh_z)) {
    absdiff_rh = abs(z_masked-rh_z[[i]])
    rh_closest_en[[names(rh_z)[[i]]]] = which(absdiff_rh==min(abs(absdiff_rh)))
  }

  # Make marks for RH based in a point
  mark = function(x, y, ...) {
    rh_arrow_default = list(x, y, x, min_z, length=.1, code = 3)
    rh_arrow_default = modifyList(rh_arrow_default, rh_arrow_args)
    do.call(arrows, rh_arrow_default)
  }

  # Find mid y for rh labels
  ymidpoint = function(x) {
    x-(x-min_z)/2
  }

  wf_default_args = list(wf, relative=FALSE, polygon=TRUE, type="l", lwd=1,
                         col="forestgreen",
                         xlab="Waveform Amplitude",
                         ylab="Elevation (m)",
                         ylim=c(ymin, ymax), xlim=qts+c(0, 0.1*abs(diff(qts))),
                         ...)
  wf_default_args = modifyList(wf_default_args, wf_args)

  do.call(plot, wf_default_args)
  par(new=TRUE)

  cumEn_default_args = list(energy_cum_normalized, z_masked, lwd=2, axes=F, bty="n", type="l", xlab = "", ylab = "", ylim=c(ymin, ymax), xlim=qts)
  cumEn_default_args = modifyList(cumEn_default_args, cumEn_args)
  do.call(plot, cumEn_default_args)

  height_axis_default = list(side=4, at = ticks, labels=ticks_label)
  height_axis_default = modifyList(height_axis_default, height_axis_args)
  do.call(axis, height_axis_default)

  height_axis_text_default = list(x=grconvertX(1, "npc"), y=grconvertY(0.5, "npc"), "Height (m)", adj=c(0.5, 5), srt=90, xpd=NA)
  height_axis_text_default = modifyList(height_axis_text_default, height_axis_text_args)
  do.call(text, height_axis_text_default)


  for (i in 2:(length(rh_z)-1)) {
    mark(energy_cum_normalized[rh_closest_en[[i]]], rh_z[[i]])
    rh_text_default = list(energy_cum_normalized[rh_closest_en[[i]]], ymidpoint(rh_z[[i]]), toupper(names(rh_z)[[i]]), pos = 2)
    rh_text_default = modifyList(rh_text_default, rh_text_args)
    do.call(text, rh_text_default)
  }

  rh100_text_default = list(qts[2]-diff(qts)/2, rh_z[[length(rh_z)]], "RH100", pos=3)
  rh100_text_default = modifyList(rh100_text_default, rh100_text_args)
  do.call(text, rh100_text_default)

  rh100_abline_default = list(rh_z[[length(rh_z)]], 0, lty="dashed")
  rh100_abline_default = modifyList(rh100_abline_default, rh100_abline_args)
  do.call(abline, rh100_abline_default)

  rh0_text_default = list(qts[2]-diff(qts)/2, rh_z[[1]], "RH0", pos=1)
  rh0_text_default = modifyList(rh0_text_default, rh0_text_args)
  do.call(text, rh0_text_default)

  rh0_abline_default = list(rh_z[[1]], 0, lty="dashed")
  rh0_abline_default = modifyList(rh0_abline_default, rh0_abline_args)
  do.call(abline, rh0_abline_default)

  gz_text_default = list(qts[1], ground_z, "GZ (m)", adj=c(0, -1))
  gz_text_default = modifyList(gz_text_default, gz_text_args)
  do.call(text, gz_text_default)

  gz_abline_default = list(ground_z, 0, lty="dashed", col="brown")
  gz_abline_default = modifyList(gz_abline_default, gz_abline_args)
  do.call(abline, gz_abline_default)


}

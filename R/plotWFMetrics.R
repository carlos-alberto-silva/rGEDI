#'GEDI-waveform plot with metrics
#'
#'@description Plots the waveform with overlaid RH metrics
#'
#'@usage plotWFMetrics(level1b, level2a, shot_number, rh=c(25, 50, 75),...)
#'
#'@param level1b A GEDI Level1B object (output of \code{\link[rGEDI:readLevel1B]{readLevel1B}} function). A S4 object of class "gedi.level1b".
#'@param level2a A GEDI Level2A object (output of \code{\link[rGEDI:readLevel2A]{readLevel2A}} function). A S4 object of class "gedi.level2a".
#'@param shot_number Shot number. A scalar represeting the shot number of a giving pulse.
#'@param rh integer vector. Specify which RH metrics to plot except rh0 and rh100, default c(25, 50, 75).
#'
#'@return Returns raster layer(s) of selected GEDI Canopy Cover and Vertical Profile Metric(s)
#'
#'@seealso https://lpdaac.usgs.gov/products/gedi02_bv001/
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
#'# specify the path to GEDI level1B and level2A data
#'level1bpath = paste0(outdir,"\\GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.h5"))
#'level2apath = paste0(outdir,"\\GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.h5"))
#'
#'# Reading GEDI level1B and level2A files
#'level1b<-readLevel1b(gedilevel1b)
#'level2a<-readLevel2a(gedilevel2a)
#'
#'shot_number = "19640521100108408"
#'
#'plotWFMetrics(level1b, level2a, shot_number, rh=c(25, 50, 75, 90))
#'}
#' 
#'@importFrom stats quantile
#'@importFrom graphics abline arrows axis mtext par text
#'@export
plotWFMetrics = function(level1b, level2a, shot_number, rh=c(25, 50, 75), ...) {
  # Extracting GEDI full-waveform for a giving shotnumber
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
  min_z = min(range_z)
  max_z = max(range_z)
  diff_z = abs(diff(range_z))
  wf_between = wf@dt[elevation %between% range_z,,]
  energy_offset = min(range_energy)
  energy_no_offset = (wf_between$rxwaveform - energy_offset)
  cumsum_energy = cumsum(rev(energy_no_offset))

  range_cumsum = range(cumsum_energy)
  range_abs_diff_cumsum = abs(diff(range_cumsum))
  energy_cum_normalized = ((cumsum_energy)/(range_abs_diff_cumsum/range_abs_diff))+energy_offset

  par(mar = c(5, 4, 4, 4) + 0.3)
  offset = diff_z*0.2
  ymin = min_z-offset
  ymax = max_z+offset
  wf_interest=wf@dt[wf@dt$elevation >= ymin & wf@dt$elevation <= ymax,]$rxwaveform
  qts=quantile(wf_interest, c(0.05, 1), type=1)


  z_masked = rev(wf_between$elevation)


  ticks = seq(min_z, max_z, length=4)
  ticks_label = format(ticks-min_z, digits = 2)

  rh_closest_en = list()
  for (i in 1:length(rh_z)) {
    absdiff_rh = abs(z_masked-rh_z[[i]])
    rh_closest_en[[names(rh_z)[[i]]]] = which(absdiff_rh==min(abs(absdiff_rh)))
  }

  # Make marks for RH based in a point
  mark = function(x, y, ...) {
    arrows(x, y, x, min_z, length=.1, code = 3)
  }

  # Find mid y for rh labels
  ymidpoint = function(x) {
    x-(x-min_z)/2
  }


  plot(wf, relative=FALSE, polygon=TRUE, type="l", lwd=2, col="forestgreen",
       xlab="Waveform Amplitude", ylab="Elevation (m)", ylim=c(ymin, ymax), xlim=qts+c(0, 0.1*abs(diff(qts))), ...)
  par(new=TRUE)
  plot(energy_cum_normalized, z_masked, lwd=2, axes=F, bty="n", type="l", xlab = "", ylab = "", ylim=c(ymin, ymax), xlim=qts)
  axis(side=4, at = ticks, labels=ticks_label)
  mtext("Height (m)", side=4, line=2)
  for (i in 2:(length(rh_z)-1)) {
    mark(energy_cum_normalized[rh_closest_en[[i]]], rh_z[[i]])
    text(energy_cum_normalized[rh_closest_en[[i]]], ymidpoint(rh_z[[i]]), toupper(names(rh_z)[[i]]), pos = 2)
  }
  text(qts[2]-diff(qts)/2, rh_z[[length(rh_z)]], "RH100", pos=3)
  abline(rh_z[[length(rh_z)]], 0, lty="dashed")
  text(qts[2]-diff(qts)/2, rh_z[[1]], "RH0", pos=1)
  abline(rh_z[[1]], 0, lty="dashed")

}

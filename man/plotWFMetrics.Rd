% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotWFMetrics.R
\name{plotWFMetrics}
\alias{plotWFMetrics}
\title{GEDI full waveform plot with metrics}
\usage{
plotWFMetrics(level1b, level2a, shot_number, rh=c(25, 50, 75),...)
}
\arguments{
\item{level1b}{A GEDI Level1B object (output of \code{\link[=readLevel1B]{readLevel1B()}} function).
An S4 object of class "gedi.level1b".}

\item{level2a}{A GEDI Level2A object (output of \code{\link[=readLevel2A]{readLevel2A()}} function).
An S4 object of class "gedi.level2a".}

\item{shot_number}{Shot number. A scalar representing the shot number of a giving pulse.}

\item{rh}{Integer vector. Specify which RH metrics to plot except rh0 and rh100, default c(25, 50, 75).}

\item{...}{Will be passed to the main plot.}
}
\value{
Nothing
}
\description{
Plots the waveform with overlaid RH metrics
}
\examples{
# specify the path to GEDI level1B and Level2A data (zip file)
outdir <- tempdir()
level1B_fp_zip <- system.file("extdata",
  "GEDI01_B_2019108080338_O01964_T05337_02_003_01_sub.zip",
  package = "rGEDI"
)

level2A_fp_zip <- system.file("extdata",
  "GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
  package = "rGEDI"
)

# Unzipping GEDI level1B data
level1Bpath <- unzip(level1B_fp_zip, exdir = outdir)
level2Apath <- unzip(level2A_fp_zip, exdir = outdir)

# Reading GEDI level1B and Level2A data (h5 file)
level1b <- readLevel1B(level1Bpath = level1Bpath)
level2a <- readLevel2A(level2Apath = level2Apath)

shot_number <- "19640521100108408"

plotWFMetrics(level1b, level2a, shot_number, rh = c(25, 50, 75, 90))

close(level1b)
close(level2a)
}
\seealso{
\url{https://lpdaac.usgs.gov/products/gedi02_bv002/}
}

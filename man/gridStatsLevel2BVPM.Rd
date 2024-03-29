% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gridStatsLevel2BVPM.R
\name{gridStatsLevel2BVPM}
\alias{gridStatsLevel2BVPM}
\title{Compute Grids with Descriptive Statistics of
GEDI derived Canopy Cover and Vertical Profile Metrics (Level2B)}
\usage{
gridStatsLevel2BVPM(level2BVPM, func, res)
}
\arguments{
\item{level2BVPM}{A GEDI Level2AM object (output of \code{\link[=getLevel2BVPM]{getLevel2BVPM()}} function).
An S4 object of class "data.table".}

\item{func}{The function(s) to be applied to each cell}

\item{res}{Spatial resolution in decimal degrees for the output stars raster layer}
}
\value{
Returns a stars raster layer(s) of selected GEDI Canopy Cover and Vertical Profile Metric(s)
}
\description{
This function computes a series of user defined descriptive statistics within
each grid cell for GEDI derived Canopy Cover and Vertical Profile Metrics (Level2B)
}
\examples{
# specify the path to GEDI level2B data (zip file)
outdir <- tempdir()
level2B_fp_zip <- system.file("extdata",
  "GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.zip",
  package = "rGEDI"
)

# Unzipping GEDI level2A data
level2Bpath <- unzip(level2B_fp_zip, exdir = outdir)

# Reading GEDI level2B data (h5 file)
level2b <- readLevel2B(level2Bpath = level2Bpath)

# Get GEDI derived Canopy Cover and Vertical Profile Metrics
level2BVPM <- getLevel2BVPM(level2b)
head(level2BVPM)

#' Define your own function
mySetOfMetrics <- function(x) {
  metrics <- list(
    min = min(x), # Min of z
    max = max(x), # Max of z
    mean = mean(x), # Mean of z
    sd = sd(x) # Sd of z
  )
  return(metrics)
}

#' Computing a serie of statistics of GEDI derived canopy cover
cover_stats <- gridStatsLevel2BVPM(
  level2BVPM = level2BVPM,
  func = mySetOfMetrics(cover),
  res = 0.005
)
plot(cover_stats)

#' Computing the max of the Total Plant Area Index only
pai_max <- gridStatsLevel2BVPM(level2BVPM = level2BVPM, func = max(pai), res = 0.005)
plot(pai_max)

#' Computing the Foliage Height Diversity Index only
fhd_mean <- gridStatsLevel2BVPM(level2BVPM = level2BVPM, func = mean(fhd_normal), res = 0.005)
plot(fhd_mean)

close(level2b)
}
\seealso{
\url{https://lpdaac.usgs.gov/products/gedi02_bv002/}
}

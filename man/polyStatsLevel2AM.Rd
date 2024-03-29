% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/polyStatsLevel2AM.R
\name{polyStatsLevel2AM}
\alias{polyStatsLevel2AM}
\title{Compute descriptive statistics of GEDI Elevation and Height Metrics}
\usage{
polyStatsLevel2AM(level2AM, func, id=NULL)
}
\arguments{
\item{level2AM}{A GEDI Level2AM object (output of \code{\link[=getLevel2AM]{getLevel2AM()}} function).
An S4 object of class "data.table".}

\item{func}{The function to be applied for computing the defined statistics}

\item{id}{A vector containing the polygon id for each GEDI observation. Default is NULL}
}
\value{
Returns an S4 object of class \link[data.table:data.table]{data.table::data.table}
Containing Statistics of GEDI level2A defined metrics
}
\description{
Computes a Series of Statistics from GEDI derived Elevation and Height Metrics (Level2A)
within a given area defined or not by a polygon
}
\examples{
# Specifying the path to GEDI level2A data (zip file)
outdir <- tempdir()
level2A_fp_zip <- system.file("extdata",
  "GEDI02_A_2019108080338_O01964_T05337_02_001_01_sub.zip",
  package = "rGEDI"
)

# Unzipping GEDI level2A data
level2Apath <- unzip(level2A_fp_zip, exdir = outdir)

# Reading GEDI level2A data (h5 file)
level2a <- readLevel2A(level2Apath = level2Apath)

# Specifying the path to shapefile
polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package = "rGEDI")

# Reading shapefile as sf object
library(sf)
polygon <- sf::st_read(polygon_filepath)

# Extracting GEDI Eleveation and Relative Metrics (level2A)
level2AM <- getLevel2AM(level2a)
head(level2AM)

# Clipping GEDI data by geometry
level2AM_clip <- clipLevel2AMGeometry(level2AM, polygon, split_by = "id")

#' Define your own function
mySetOfMetrics <- function(x) {
  metrics <- list(
    min = min(x), # Min of x
    max = max(x), # Max of x
    mean = mean(x), # Mean of x
    sd = sd(x) # Sd of x
  )
  return(metrics)
}

# Computing the maximum of RH100
RH100max <- polyStatsLevel2AM(level2AM_clip, func = max(rh100), id = NULL)

# Computing the maximum of RH100 stratified by polygon
RH100max_poly <- polyStatsLevel2AM(level2AM_clip, func = max(rh100), id = NULL)

# Computing a serie statistics for GEDI metrics stratified by polygon
RH100metrics <- polyStatsLevel2AM(level2AM_clip,
  func = mySetOfMetrics(rh100),
  id = level2AM_clip$id
)

head(RH100metrics)

close(level2a)
}
\seealso{
\url{https://lpdaac.usgs.gov/products/gedi02_av002/}
}

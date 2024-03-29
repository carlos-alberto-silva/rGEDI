% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clipLevel2BPAIProfile.R
\name{clipLevel2BPAIProfileGeometry}
\alias{clipLevel2BPAIProfileGeometry}
\title{Clip GEDI Plant Area Index profile by geometry}
\usage{
clipLevel2BPAIProfileGeometry(level2BPAIProfile, polygon, split_by = NULL)
}
\arguments{
\item{level2BPAIProfile}{A GEDI Level2B object (output of \code{\link[=getLevel2BPAIProfile]{getLevel2BPAIProfile()}} function).
An S4 object of class "data.table".}

\item{polygon}{Polygon. An object of class \code{\link[sf:sf]{sf::sf}},
which can be loaded as an ESRI shapefile using \link[sf:st_read]{sf::st_read} function in the \emph{sf} package.}

\item{split_by}{Polygon id. If defined, GEDI data will be clipped by each polygon using the
attribute specified by \code{split_by} from the attribute table.}
}
\value{
Returns an S4 object of class \link[data.table:data.table]{data.table::data.table}
containing the Plant Area Index profile data.
}
\description{
This function clips GEDI level2B derived
Plant Area Index profile within a given geometry
}
\examples{
# Specifying the path to GEDI level2B data (zip file)
outdir <- tempdir()
level2B_fp_zip <- system.file("extdata",
  "GEDI02_B_2019108080338_O01964_T05337_02_001_01_sub.zip",
  package = "rGEDI"
)

# Unzipping GEDI level2A data
level2Bpath <- unzip(level2B_fp_zip, exdir = outdir)

# Reading GEDI level2B data (h5 file)
level2b <- readLevel2B(level2Bpath = level2Bpath)

# Extracting GEDI Plant Area Index profile
level2BPAIProfile <- getLevel2BPAIProfile(level2b)

# Specifying the path to shapefile
polygon_filepath <- system.file("extdata", "stands_cerrado.shp", package = "rGEDI")

# Reading shapefile as sf object
library(sf)
polygon <- sf::st_read(polygon_filepath)

# Clipping GEDI Plant Area Index profile by geometry
level2b_clip_geometry <- clipLevel2BPAIProfileGeometry(
  level2BPAIProfile,
  polygon,
  split_by = "id"
)

close(level2b)
}
\seealso{
\url{https://lpdaac.usgs.gov/products/gedi02_bv002/}
}

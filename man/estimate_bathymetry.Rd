% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_bathymetry.R
\name{estimate_bathymetry}
\alias{estimate_bathymetry}
\title{Estimate bathymetry from a shoreline polygon using distance transformation.}
\usage{
estimate_bathymetry(shoreline, max_depth, res = 10)
}
\arguments{
\item{shoreline}{sf object of lake shoreline.}

\item{max_depth}{numeric. The maximum depth of the lake.}

\item{res}{numeric resolution of output raster in metres.}
}
\value{
SpatRaster object with the estimated bathymetry.
}
\description{
This function estimates the bathymetry of a lake from a shoreline polygon using
a distance transformation. The function calculates the distance from the
shoreline polygon to the edge of the lake and uses this distance to estimate
the bathymetry. The maximum depth of the lake is provided as an input parameter.
This is the formula used in the Hollister and Milstead (2010) paper.
}
\examples{
shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
package = "bathytools"))
max_depth <- 81.6
bathy <- estimate_bathymetry(shoreline = shoreline, max_depth = max_depth)
terra::plot(bathy)
}
\references{
Hollister, J. W., W.B. Milstead (2010). Using GIS to Estimate
Lake Volume from Limited Data. Lake and Reservoir Management.
26(3)194-199.
\doi{10.1080/07438141.2010.504321}
}

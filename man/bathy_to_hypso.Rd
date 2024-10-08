% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bathy_to_hypso.R
\name{bathy_to_hypso}
\alias{bathy_to_hypso}
\title{Convert bathymetry raster to hypsograph}
\usage{
bathy_to_hypso(bathy_raster, surface = 0, depths = 1)
}
\arguments{
\item{bathy_raster}{SpatRaster object with the bathymetry data.}

\item{surface}{numeric. The surface elevation of the lake. Default is 0.}

\item{depths}{numeric. The depths at which to calculate the area. If a single
numeric value is provided, the function will calculate the area at each depth
from the surface to the minimum depth of the bathymetry raster at intervals of
the provided value. If a vector of numeric values is provided, the function
will calculate the area at each depth specified in the vector. Default is 1.}
}
\value{
data.frame
}
\description{
This function converts a bathymetry raster to a hypsograph,
which is a data.frame with the area of the lake at different depths. The
function takes a bathymetry raster as input and returns a data frame with the
depth and area at each depth.
}
\examples{
shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
package = "bathytools"))
point_data <- readRDS(system.file("extdata/depth_points.rds",
package = "bathytools"))
bathy_raster <- rasterise_bathy(shoreline = shoreline,
point_data = point_data, crs = 2193)
hyps <- bathy_to_hypso(bathy_raster = bathy_raster)
}

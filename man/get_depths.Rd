% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_depths.R
\name{get_depths}
\alias{get_depths}
\title{Get depths for hypsograph or contours}
\usage{
get_depths(bathy_raster, surface = 0, depths)
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
numeric
}
\description{
This function calculates the depths at which to calculate the
area for a hypsograph or contours. If a single numeric value is provided, the
function will calculate the area at each depth from the surface to the minimum
depth of the bathymetry raster at intervals of the provided value. If a vector
of numeric values is provided, the function will calculate the area at each
depth specified in the vector. If the minimum depth of the bathymetry raster
is not in the vector of depths, it will be added to the vector.
}

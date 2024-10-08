#' @title Get depths for hypsograph or contours
#' @description This function calculates the depths at which to calculate the
#' area for a hypsograph or contours. If a single numeric value is provided, the
#' function will calculate the area at each depth from the surface to the minimum
#' depth of the bathymetry raster at intervals of the provided value. If a vector
#' of numeric values is provided, the function will calculate the area at each
#' depth specified in the vector. If the minimum depth of the bathymetry raster
#' is not in the vector of depths, it will be added to the vector.
#' @param bathy_raster SpatRaster object with the bathymetry data.
#' @param surface numeric. The surface elevation of the lake. Default is 0.
#' @param depths numeric. The depths at which to calculate the area. If a single
#' numeric value is provided, the function will calculate the area at each depth
#' from the surface to the minimum depth of the bathymetry raster at intervals of
#' the provided value. If a vector of numeric values is provided, the function
#' will calculate the area at each depth specified in the vector. Default is 1.
#' @return numeric
#' @export
#' @importFrom terra minmax
get_depths <- function(bathy_raster, surface = 0, depths) {

  if (!is(bathy_raster, "SpatRaster")) {
    stop("bathy_raster must be a SpatRaster object")
  }

  # Get min/max of bathy object
  mm <- terra::minmax(bathy_raster)

  if (length(depths) == 1) {
    depth_out <- seq(from = surface, to = mm[1], by = -depths)
  } else if (length(depths) > 1) {
    depth_out <- depths
  }
  if (!mm[1] %in% depths) {
    depth_out <- c(depth_out, mm[1])
  }

  # Round depth_out to 2 decimal
  depth_out <- round(depth_out, 2)
  return(depth_out)
}

#' Get Lake Depth from Bathymetry Raster
#' This function calculates the lake depth from a bathymetric raster.
#' @param bathy_raster SpatRaster object with the bathymetry data.
#' @return numeric. The lake depth.
#' @export
#' @importFrom terra minmax
#'

get_lake_depth <- function(bathy_raster) {
  # Check if the input is a SpatRaster object
  if (!inherits(bathy_raster, "SpatRaster")) {
    stop("Input must be a SpatRaster object.")
  }
  mm <- terra::minmax(bathy_raster)
  lake_depth <- abs(mm[1])
  return(lake_depth)
}

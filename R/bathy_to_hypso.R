#' @title Convert bathymetry raster to hypsograph
#'
#' @description This function converts a bathymetry raster to a hypsograph,
#' which is a data.frame with the area of the lake at different depths. The
#' function takes a bathymetry raster as input and returns a data frame with the
#' depth and area at each depth.
#'
#' @inheritParams merge_bathy_dem
#' @param surface numeric. The surface elevation of the lake. Default is 0.
#' @param depths numeric. The depths at which to calculate the area. If a single
#' numeric value is provided, the function will calculate the area at each depth
#' from the surface to the minimum depth of the bathymetry raster at intervals of
#' the provided value. If a vector of numeric values is provided, the function
#' will calculate the area at each depth specified in the vector. Default is 1.
#'
#' @importFrom terra minmax res values
#' @importFrom dplyr arrange
#'
#' @return data.frame
#' @export
#'
#' @examples
#' shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
#' package = "bathytools"))
#' point_data <- readRDS(system.file("extdata/depth_points.rds",
#' package = "bathytools"))
#' bathy_raster <- rasterise_bathy(shoreline = shoreline,
#' point_data = point_data, crs = 2193)
#' hyps <- bathy_to_hypso(bathy_raster = bathy_raster)

bathy_to_hypso <- function(bathy_raster, surface = 0, depths = 1) {

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

  # Get resolution of bathy object
  res <- terra::res(bathy_raster)

  # Calculate areas of each depth
  areas <- rep(NA, length(depth_out))
  areas <- sapply(depth_out, \(d) {
    sum(terra::values(bathy_raster) <= d, na.rm = TRUE) * res[1] * res[2]
  })
  df <- data.frame(depth = depth_out, area = areas) |>
    dplyr::arrange(desc(depth))
  return(df)
}

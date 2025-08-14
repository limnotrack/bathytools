#' @title Convert bathymetry raster to hypsograph
#'
#' @description This function converts a bathymetry raster to a hypsograph,
#' which is a data.frame with the area of the lake at different depths. The
#' function takes a bathymetry raster as input and returns a data frame with the
#' depth and area at each depth.
#'
#' @inheritParams merge_bathy_dem
#' @inheritParams get_depths
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

bathy_to_hypso <- function(bathy_raster, surface = 0, depths = NULL) {

  # Get depths
  depth_out <- get_depths(bathy_raster, surface, depths)


  # Get resolution of bathy object
  res <- terra::res(bathy_raster)

  # Calculate areas of each depth
  areas <- rep(NA, length(depth_out))
  areas <- sapply(depth_out, \(d) {
    sum(terra::values(bathy_raster) <= d, na.rm = TRUE) * res[1] * res[2]
  })
  df <- data.frame(depth = depth_out, area = areas) |>
    dplyr::mutate(elev = depth) |> 
    dplyr::arrange(desc(depth)) |> 
    dplyr::select(elev, depth, area)
  return(df)
}

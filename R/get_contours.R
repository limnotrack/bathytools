#' Get contours from a bathymetric raster
#'
#' @inheritParams bathy_to_hypso
#'
#' @return sf LINESTRING of the contours.
#' @export
#'
#' @importFrom terra as.contour values
#' @importFrom dplyr rename
#' @importFrom sf st_as_sf
#'
#' @examples
#' shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
#' package = "bathytools"))
#' point_data <- readRDS(system.file("extdata/depth_points.rds",
#' package = "bathytools"))
#' bathy_raster <- rasterise_bathy(shoreline = shoreline,
#' point_data = point_data, crs = 2193, res = 8)
#' contours <- get_contours(bathy_raster = bathy_raster)
get_contours <- function(bathy_raster, surface = 0, depths = 1) {

  # Get depths
  depth_out <- get_depths(bathy_raster, surface, depths)

  # Cast to linestring
  # shoreline <- sf::st_cast(shoreline$geometry, "LINESTRING") |>
  #   sf::st_sf(depth = 0)
  # sf::st_geometry(shoreline) <- "geometry"
  # Get contours
  cont <- terra::as.contour(x = bathy_raster, levels = depth_out) |>
    sf::st_as_sf() |>
    dplyr::rename(depth = level)

  # Check if min depth is in depths
  if (!min(terra::values(bathy_raster)) %in% depth_out) {
    warning("Minimum depth not in depths.")
    # cont <- dplyr::bind_rows(cont, tibble::tibble(depth = min(terra::values(bathy_raster))))
  }

  return(cont)

}

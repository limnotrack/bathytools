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
get_contours <- function(bathy_raster, surface = 0, depths = NULL) {

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
  mm <- terra::minmax(bathy_raster)
  if (!mm[1] %in% depth_out) {
    warning("Minimum depth not in depths. Adding minimum depth to contours.")

    min_contour <- min(cont$depth)

    btm <- bathy_raster <= (mm[1] + 0.1)
    # terra::plot(btm)
    btm_cont <- btm |>
      terra::as.polygons() |>
      sf::st_as_sf() |>
      dplyr::filter(depth == 1) |>
      dplyr::mutate(depth = mm[1])

    cont <- dplyr::bind_rows(cont, btm_cont)
  }

  return(cont)

}

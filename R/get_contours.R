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
#' depth_points <- readRDS(system.file("extdata/depth_points.rds",
#' package = "bathytools"))
#' bathy_raster <- rasterise_bathy(shoreline = shoreline,
#' depth_points = depth_points, crs = 2193, res = 8)
#' contours <- get_contours(bathy_raster = bathy_raster)
get_contours <- function(bathy_raster, surface = 0, depths = NULL) {

  # Get depths
  depth_out <- get_depths(bathy_raster, surface, depths)

  mm <- terra::minmax(bathy_raster)
  # If minimum depth is less than 0, then we need to adjust the output depths to
  # negative
  if (mm[2] < 0) {
    depth_out <- -depth_out
  }
  cont <- terra::as.contour(x = bathy_raster, levels = depth_out) |>
    sf::st_as_sf() |>
    dplyr::rename(depth = level)

  # Check if min depth is in depths
  mm <- terra::minmax(bathy_raster) |> 
    round(2)
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

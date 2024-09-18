#' Generate depth points for interpolation
#'
#' @inheritParams generate_depth_points
#' @inheritParams interpolate_points
#' @param method character interpolation method. Options are 'MBA' (default),
#' @inheritParams sf::st_transform
#' @inheritParams MBA::mba.surf
#'
#' @return sf object of depth points.
#' @export
#'
#' @examples
#' shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
#' package = "bathytools"))
#' point_data <- readRDS(system.file("extdata/depth_points.rds",
#' package = "bathytools"))
#' bathy <- rasterise_bathy(shoreline = shoreline, point_data = point_data,
#' crs = 2193)
#'

rasterise_bathy <- function(shoreline, islands = NULL, point_data = NULL,
                            contours = NULL, res = 2, subsample = TRUE,
                            crs, method = "MBA", print_plot = TRUE,
                            n = 1, m = 1, h = 8) {

  all <- generate_depth_points(shoreline = shoreline, islands = islands,
                               point_data = point_data, contours = contours,
                               res = res, subsample = subsample, crs = crs)

  bathy <- interpolate_points(point_data = all, shoreline = shoreline,
                              res = res, method = method, n = n, m = m, h = h,
                              crs = crs, print_plot = print_plot)
  return(bathy)
}




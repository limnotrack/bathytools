#' Generate depth points for interpolation
#'
#' @inheritParams generate_depth_points
#' @inheritParams interpolate_points
#' @param method character interpolation method. Options are nearest neighbour 
#' ('nn') (default), multilevel B-splines ('MBA'), thin plate spline ('tps'),
#' and inverse distance weighting ('idw').
#' @inheritParams sf::st_transform
#' @inheritParams MBA::mba.surf
#'
#' @return sf object of depth points.
#' @export
#'
#' @examples
#' shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
#' package = "bathytools"))
#' depth_points <- readRDS(system.file("extdata/depth_points.rds",
#' package = "bathytools"))
#' bathy <- rasterise_bathy(shoreline = shoreline, depth_points = depth_points,
#' crs = 2193)
#'

rasterise_bathy <- function(shoreline, islands = NULL, depth_points = NULL,
                            contours = NULL, res = 2, subsample = TRUE,
                            crs, method = "nn", print_plot = TRUE,
                            n = 1, m = 1, h = 8) {

  depth_points <- generate_depth_points(shoreline = shoreline, 
                                        islands = islands,
                                        depth_points = depth_points, 
                                        contours = contours, res = res,
                                        subsample = subsample, crs = crs)

  bathy <- interpolate_points(depth_points = depth_points, 
                              shoreline = shoreline, islands = islands,
                              res = res, method = method, n = n, m = m, h = h,
                              crs = crs, print_plot = print_plot)
  return(bathy)
}




#' Generate depth points for interpolation
#'
#' @inheritParams generate_depth_points
#' @inheritParams interpolate_points
#' @param method character interpolation method. Options are c('MBA', 'tps',
#' 'nn', 'idw'). Default is 'nn' (nearest neighbor). 'MBA' uses the
#' MBA::mba.surf function, 'tps' uses a thin plate spline from the fields
#' package, 'nn' uses nearest neighbor interpolation from the terra package,
#' and 'idw' uses inverse distance weighting from the terra package.
#' @inheritParams sf::st_transform
#' @inheritParams MBA::mba.surf
#' 
#' @importFrom cli cli_progress_step
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
                            crs, method = c("MBA", "tps", "nn", "idw"),
                            print_plot = TRUE,
                            n = 1, m = 1, h = 8) {

  method <- rlang::arg_match(method)
  
  missing_islands <- detect_islands(shoreline)
  if (!is.null(missing_islands)) {
    if (is.null(islands)) {
      islands <- missing_islands
    } else {
      islands <- dplyr::bind_rows(islands, missing_islands)
    }
  }
  shoreline <- detect_shoreline(shoreline)
  
  # cli::cli_progress_step("Generating depth points for interpolation")
  depth_points <- generate_depth_points(shoreline = shoreline, 
                                        islands = islands,
                                        depth_points = depth_points, 
                                        contours = contours, res = res,
                                        subsample = subsample, crs = crs)

  # cli::cli_progress_step("Interpolating depth points to raster")
  bathy <- interpolate_points(depth_points = depth_points, 
                              shoreline = shoreline, 
                              islands = islands,
                              res = res, method = method, n = n, m = m, h = h,
                              crs = crs, print_plot = print_plot)
  if (!is.null(islands)) {
    bathy <- terra::mask(bathy, islands, inverse = TRUE)
  }
  return(bathy)
}




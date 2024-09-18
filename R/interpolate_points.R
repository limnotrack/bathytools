#' Interpolate point depth data to a raster
#'
#' @inheritParams generate_depth_points
#' @param method character interpolation method. Options are 'MBA' (default).
#' @inheritParams MBA::mba.surf
#' @param print_plot logical print plot of interpolated raster.
#'
#' @importFrom MBA mba.surf
#' @importFrom terra rast resample crs interpolate plot mask res
#' @importFrom sf st_transform st_crs st_make_grid st_coordinates st_bbox
#' @importFrom dplyr rename mutate filter case_when
#'
#' @return SpatRaster
#' @export
#'
#' @examples
#' shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
#' package = "bathytools"))
#' point_data <- readRDS(system.file("extdata/depth_points.rds",
#'                                   package = "bathytools"))
#' depth_points <- generate_depth_points(shoreline = shoreline,
#' point_data = point_data)
#' bathy <- interpolate_points(point_data = depth_points, shoreline = shoreline,
#' crs = 2193)
#'

interpolate_points <- function(point_data, shoreline, crs, res = 2,
                               method = "MBA", n = 1, m = 1, h = 8,
                               print_plot = TRUE) {

  message("Interpolating to raster... [", format(Sys.time()), "]")

  # Convert points to a data.frame
  coords <- point_data |>
    sf::st_transform(crs) |>
    # Add x and y columns
    sf::st_coordinates() |>
    as.data.frame() |>
    dplyr::rename(x = X, y = Y) |>
    dplyr::mutate(z = point_data$depth)


  # Transform point data to sf object and rename depth to z
  data_sf <- point_data |>
    sf::st_transform(crs) |>
    dplyr::rename(z = depth)

  # Store CRS for later
  point_crs <- terra::crs(data_sf)
  bbox <- sf::st_bbox(data_sf)

  # If shoreline CRS is different from point data, transform
  if (sf::st_crs(shoreline) != crs) {
    shoreline <- sf::st_transform(shoreline, crs)
  }

  # Create a buffer around the lake and generate a grid
  buff_lake <- sf::st_buffer(shoreline, 5 * res)
  grd_sf <- sf::st_make_grid(x = buff_lake,
                             cellsize = res,
                             what = "centers", square = TRUE)

  # Convert the grid into an xyz dataframe
  grd <- grd_sf |>
    sf::st_coordinates() |>
    as.data.frame() |>
    dplyr::mutate(z = 0)

  # Make a raster grid
  grd_ras <- grd |>
    terra::rast(crs = point_crs)

  if (method == "MBA") {

    mba <- MBA::mba.surf(coords, no.X = (bbox[["xmax"]]-bbox[["xmin"]])/res +1,
                         no.Y = (bbox[["xmax"]]-bbox[["xmin"]])/res+1,
                         n = n, m = m, h = h)
    # Create a data frame of all combinations of x and y
    df_mba <- expand.grid(x = mba$xyz.est$x, y = mba$xyz.est$y)

    # Flatten the z values and add them to the data frame
    df_mba$z <- as.vector(mba$xyz.est$z)
    df_mba <- df_mba |>
      dplyr::filter(!is.na(z))

    # Convert dataframe to raster
    interp <- terra::rast(df_mba, crs = point_crs)
    # Resample to ensure it is the correct resolution
    interp2 <- interp
    terra::res(interp2) <- res
    resample <- terra::resample(interp, interp2)
    interp <- resample

  }
  # Mask the raster with the shoreline
  bathy <- terra::mask(interp, shoreline)

  # Check max values and ensure all are below 0
  mm <- terra::minmax(bathy)
  if (mm[2] >= 0) {
    message("Adjusting depths >= 0")
    # surf <- bathy > -0.5
    # terra::plot(surf)
    bathy[bathy > -0.5] <- -0.5
  }

  message("Finished! [", format(Sys.time()), "]")


  # Rename the variable in the raster
  names(bathy) <- "depth"

  if (print_plot) {
    terra::plot(bathy)
  }
  return(bathy)
}



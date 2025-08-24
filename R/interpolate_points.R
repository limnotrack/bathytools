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
#' package = "bathy_rastertools"))
#' depth_points <- readRDS(system.file("extdata/depth_points.rds",
#'                                   package = "bathy_rastertools"))
#' depth_points <- generate_depth_points(shoreline = shoreline,
#' depth_points = depth_points)
#' bathy_raster <- interpolate_points(depth_points = depth_points, shoreline = shoreline,
#' crs = 2193)
#'

interpolate_points <- function(depth_points, shoreline, islands = NULL, crs,
                               res = 2, method = "nn", n = 1, m = 1, h = 8,
                               print_plot = TRUE) {

  message("Interpolating to raster... [", format(Sys.time()), "]")

  # Convert points to a data.frame
  coords <- depth_points |>
    sf::st_transform(crs) |>
    # Add x and y columns
    sf::st_coordinates() |>
    as.data.frame() |>
    dplyr::rename(x = X, y = Y) |>
    dplyr::mutate(z = depth_points$depth)


  # Transform point data to sf object and rename depth to z
  data_sf <- depth_points |>
    sf::st_transform(crs) |>
    dplyr::rename(z = depth)
  
  # Convert to a vector for interpolation
  data_vec <- terra::vect(data_sf)

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

    mba <- MBA::mba.surf(coords, 
                         no.X = (bbox[["xmax"]]-bbox[["xmin"]]) / (res + 1),
                         no.Y = (bbox[["xmax"]]-bbox[["xmin"]]) / (res + 1),
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
    interp <- resample |> 
      terra::mask(shoreline)
    # terra::plot(interp)

  } else if (method == "tps") {
    
    # Wrapper to make sure terra passes the right format
    tps_fun <- function(model, newdata) {
      # newdata from terra has cols x and y
      new_coords <- as.matrix(newdata[, c("x", "y")])
      predict(model, new_coords)
    }
    
    # Thin plate spline model
    # idx <- sample(1:nrow(coords), size = (0.75), replace = FALSE)
    idx <- 1:nrow(coords)
    tps <- fields::Tps(x = coords[idx, 1:2], Y = coords$z[idx])
    tps <- fields::fastTps(x = coords[, 1:2], Y = coords$z, aRange = 2)
    
    
    # Interpolate to raster
    interp <- terra::interpolate(grd_ras, model = tps, fun = tps_fun) 
    terra::plot(interp, main = "Thin Plate Spline Interpolation")
    

    # tps <- fields::fastTps(x = coords[, 1:2], Y = coords$z, aRange = 2)
    # idx <- sample(1:nrow(coords), size = 500, replace = FALSE)
    tps <- fields::Tps(x = coords[idx, 1:2], Y = coords$z[idx])
    interp <- terra::interpolate(object = grd_ras, model = tps)
    
    # gstat::krige(z ~ 1, locations = data_sf, newdata = grd_ras,
    #              nmax = 12, set = list(idp = 2))
    
    i <- idw(z~1, data_sf, grd)
    # Set up gstat model
    idw_model <- gstat::gstat(formula = z ~ 1,
                              locations = data_sf,
                              nmax = 12,          # number of neighbors
                              set = list(idp = 2))  # inverse distance power
    
    # Predict to grid
    interp <- terra::interpolate(grd_ras, idw_model)
    
  } else if (method == "nn") {
    
    # Convert raster cell centers to points
    grid_pts <- terra::as.points(grd_ras)
    
    # Nearest distance from each grid point to your survey points
    dists <- terra::distance(grd_ras, data_vec)
    radius <- max(terra::values(dists), na.rm = TRUE)
    
    interp <- terra::interpNear(x = grd_ras, y = data_vec, field = "z",
                               radius = radius, 
                               interpolate = TRUE) 
    # terra::plot(interp, main = "Nearest Neighbor Interpolation")
    
    
  } else if (method == "idw") {
    interp <- terra::interpIDW(x = grd_ras, y = data_vec, field = "z",
                               minPoints = 3, maxPoints = 10, smooth = 10,
                               radius = 1400) |> 
      terra::mask(shoreline)
  }
  # Mask the raster with the shoreline
  bathy_raster <- terra::mask(interp, shoreline)
  if (!is.null(islands)) {
    islands <- sf::st_transform(islands, crs)
    bathy_raster <- terra::mask(bathy_raster, islands, inverse = TRUE)
  }

  # Check max values and ensure all are below 0
  mm <- terra::minmax(bathy_raster)
  if (mm[2] >= 0) {
    min_depth <- mm[1]
    adj_depth <- round(0.01 * min_depth, 2)
    adj_depth <- ifelse(adj_depth > -0.4, -0.4, adj_depth)
    message(paste("Adjusting depths >= 0 to ", adj_depth, "m"))

    # surf <- bathy_raster > -0.5
    # terra::plot(surf)
    bathy_raster[bathy_raster > adj_depth] <- adj_depth
  }

  message("Finished! [", format(Sys.time()), "]")


  # Rename the variable in the raster
  names(bathy_raster) <- "depth"

  if (print_plot) {
    # depths <- unique(depth_points$depth)
    # hyps <- bathytools::bathy_to_hypso(bathy_raster, depths = depths)
    # plot(hyps$area, hyps$depth)
    terra::plot(bathy_raster)
    
  }
  return(bathy_raster)
}



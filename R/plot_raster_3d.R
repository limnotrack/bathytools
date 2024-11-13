#' Plot a 3D raster
#'
#' @param x SpatRaster object
#' @inheritParams terra::aggregate
#' @param mask_to_lake logical, mask the raster to the lake. Default is FALSE.
#' @inheritParams generate_depth_points
#' @inheritParams tm_dem_bath
#' @param aspectratio_z numeric, aspect ratio of the z axis. Default is 0.2.
#' @param dim_px numeric, dimension of the plot in pixels. Default is 300.
#'
#' @importFrom terra ext mask aggregate as.matrix values
#' @importFrom sf st_area st_buffer
#' @importFrom plotly plot_ly layout add_surface add_trace
#' @importFrom units drop_units
#' @importFrom dplyr rename mutate
#'
#' @return plotly object
#' @export
#'
#' @examples
#' r <- terra::rast(system.file("ex/elev.tif", package="terra"))
#' plot_raster_3d(r)

plot_raster_3d <- function(x, fact = 0, cols = "topo", split_lake = FALSE,
                           mask_to_lake = FALSE,
                           shoreline = NULL, lake_elev = NULL, add_shoreline = FALSE,
                           aspectratio_z = 0.2, dim_px = 300) {

  # Colours for the terrain
  terr_cols <- c("#C7E9C0", "#41AB5D", "#E6E600", "#EBB25E", "#F0C9C0")
  bathy_cols <- c("#084594", "white")
  if (cols == "topo") {
    base_cols <- terr_cols
    col_label = "Topography"
  } else if (cols == "bathy") {
    base_cols <- bathy_cols
    col_label = "Bathymetry"
  }


  # Check factor for aggregation
  if (fact != 0) {
    if (fact < 0) {
      stop("fact must be greater than 0")
    }
  }

  # Get CRS from the raster
  crs <- sf::st_crs(x)
  # Set xaxis and yaxis labels
  if (crs$input == "WGS 84") {
    xlab <- "Longitude"
    ylab <- "Latitude"
  } else {
    xlab <- "Easting"
    ylab <- "Northing"
  }

  # Set up the plotly scene ----
  scene <- list(camera = list(eye = list(x = -0.5, y = -0.75, z = 1.25)),
               aspectmode = "manual", aspectratio = list(x = 1, y = 1,
                                                         z = aspectratio_z),
               xaxis = list(title = xlab),
               yaxis = list(title = ylab),
               zaxis = list(title = "Elevation (masl)", nticks = 4))


  p <- plotly::plot_ly() |>
    plotly::layout(scene = scene)

  # make the surfaces roughly the same size (per 'dim' arg)
  if (fact != 0) {
    ras_agg <- terra::aggregate(x, fact = fact)
  } else {
    ras_agg <- x
  }
  reso <- terra::res(ras_agg)
  ext_agg <- terra::ext(ras_agg)
  buff <- diff(ext_agg[1:2]) / dim_px * 1.5

  if (!split_lake) {
    # Calculate area of dem
    # ext_dem  <- terra::ext(x)
    # area_dem <- diff(ext_dem[1:2]) * diff(ext_dem[3:4])
    # extent of the aggregation

    # matrix for display
    ras_mat <- ras_agg |>
      terra::as.matrix(x = _, wide = TRUE)
    # reverse the matrix in y dimension
    # ras_mat <- ras_mat[, ncol(ras_mat):1]

    # add it to the surface
    p <- p |>
      plotly::add_surface(x = seq(ext_agg[1], ext_agg[2],
                                  length.out = ncol(ras_mat)),
                          y = seq(ext_agg[4], ext_agg[3],
                                  length.out = nrow(ras_mat)),
                          z = ~ras_mat,
                          colors = terr_cols,
                          colorbar = list(title = col_label),
                          colorscale = list(seq(0, 1,
                                                length.out = length(base_cols)),
                                            base_cols))

  } else {
    # mask the lake off
    if (is.null(shoreline)) {
      stop("shoreline must be provided if using split_lake = TRUE")
    }
    bathy_agg <- terra::mask(ras_agg, shoreline, touches = FALSE)
    if (!is.null(lake_elev)) {
      bathy_agg[bathy_agg > lake_elev] <- NA
    }
    # create a false lip around the lake
    bathy_agg[is.na(bathy_agg[])] <- max(terra::values(bathy_agg), na.rm = TRUE)
    bathy_agg <- terra::mask(bathy_agg, sf::st_buffer(shoreline, buff * 2),
                             touches = FALSE)
    bathy_mat <- bathy_agg |>
      terra::as.matrix(x = _, wide = TRUE)


    # process the catchment if needed
    if (!is.null(shoreline)) {
      land_mat <- ras_agg |>
        terra::mask(sf::st_buffer(shoreline, -buff), touches = FALSE,
                    inverse = TRUE) |>
        terra::as.matrix(x = _, wide = T)
      # land_mat <- land_mat[, ncol(land_mat):1]
    }
    if (!mask_to_lake) {
      # add it to the surface
      p <- p |>
        plotly::add_surface(x = seq(ext_agg[1], ext_agg[2],
                                    length.out = ncol(land_mat)),
                            y = seq(ext_agg[4], ext_agg[3],
                                    length.out = nrow(land_mat)),
                            z = ~land_mat,
                            colors = terr_cols,
                            colorbar = list(title = "Topography"),
                            colorscale =
                              list(seq(0, 1, length.out = length(terr_cols)),
                                   terr_cols))
    }
    ext_bathy_agg <- terra::ext(bathy_agg)
    p <- p |>
      plotly::add_surface(x = seq(ext_bathy_agg[1], ext_bathy_agg[2],
                                  length.out = ncol(bathy_mat)),
                          y = seq(ext_bathy_agg[4], ext_bathy_agg[3],
                                  length.out = nrow(bathy_mat)),
                          z = ~bathy_mat,
                          colors = bathy_cols,
                          colorbar = list(title = "Bathymetry"),
                          colorscale = list(c(0, 1), bathy_cols))

  }

  # Convert shoreline to trace
  if (!is.null(shoreline) & add_shoreline) {
    if (is.null(lake_elev)) {
      lake_elev <- 0
    }

    shore_trace <- sf::st_coordinates(shoreline) |>
      as.data.frame() |>
      dplyr::rename(x = X, y = Y) |>
      dplyr::mutate(elevation = lake_elev)
    p <- p |>
      plotly::add_trace(x = shore_trace$x, y = shore_trace$y,
                        z = shore_trace$elevation,
                        type = "scatter3d", mode = "lines",
                        line = list(color = "orange", width = reso[1]))
  }




  # p
  return(p)

}


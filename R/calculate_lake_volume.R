#' Calculate the volume of a lake using bathymetry data or a hypsograph
#'
#' @param bathy_raster SpatRaster object with the bathymetry data.
#' @param hyps data.frame with columns 'depth' and 'area'.
#' @param depth numeric. The depth to which to calculate the volume. If
#' provided, the volume will be calculated to this depth. If not provided, the
#' volume will be calculated to the maximum depth of the bathymetry raster or
#' the hypsograph.
#' @param return_rast logical. If TRUE, return a raster with the calculated
#' volume in each grid cell. Default is FALSE.
#'
#' @return numeric. The volume of the lake in cubic meters.
#' @export
#' @examples
#' shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
#' package = "bathytools"))
#' point_data <- readRDS(system.file("extdata/depth_points.rds",
#' package = "bathytools"))
#' bathy_raster <- rasterise_bathy(shoreline = shoreline,
#' point_data = point_data, crs = 2193)
#' calculate_lake_volume(bathy_raster = bathy_raster)
#'

calculate_lake_volume <- function(bathy_raster = NULL, hyps = NULL,
                                  depth = NULL, return_rast = FALSE) {

  if (!is.null(bathy_raster)) {
    if (!is(bathy_raster, "SpatRaster")) {
      stop("Input must be a SpatRaster object")
    }

    mm <- terra::minmax(bathy_raster)
    if (!is.null(depth)) {
      upper <- mm[1] + depth
    } else {
      upper <- 0
    }

    adj_bathy <- terra::clamp(bathy_raster, upper = upper, values = FALSE)
    # Calculate the area of each raster cell (in square meters)
    cell_area <- terra::cellSize(adj_bathy, unit = "m")

    # Calculate volume for each cell (cell area * depth)
    # Make sure to take the absolute value if depth is negative (e.g., below sea level)
    cell_volume <- cell_area * abs(terra::values(adj_bathy))

    if (return_rast) {
      return(cell_volume)
    }

    # Sum the cell volumes to get the total lake volume (in cubic meters)
    total_volume <- sum(terra::values(cell_volume), na.rm = TRUE)
  } else if (!is.null(hyps)) {
    if (!is.data.frame(hyps)) {
      stop("Input must be a data.frame object")
    }
    if (!all(names(hyps) %in% c("depth", "area"))) {
      stop("Input data.frame must have columns 'depth' and 'area'")
    }
    if (!is.null(depth)) {
      upper <- min(hyps$depth) + depth
    } else {
      upper <- 0
    }
    # Fit a linear model
    model <- lm(area ~ depth, data = hyps)

    hyps <- hyps |>
      dplyr::filter(depth <= upper) |>
      dplyr::arrange(dplyr::desc(depth))

    if (!upper %in% hyps$depth) {
      # Predict new values, including extrapolated ones
      new_x <- data.frame(depth = c(upper))
      surf_area <- predict(model, newdata = new_x)
      hyps <- dplyr::bind_rows(data.frame(depth = upper, area = surf_area), hyps)
    }

    # Calculate the incremental volumes
    hyps$volume <- with(hyps, {
      # Calculate volumes between successive depth intervals
      vol <- (area[-nrow(hyps)] + area[-1]) / 2 * diff(abs(depth))
      c(vol, 0)  # Add a zero for the last depth as it has no next depth
    })

    # Sum the incremental volumes to get total lake volume
    total_volume <- sum(hyps$volume)
  } else {
    stop("Either bathy_raster or hyps must be provided")
  }

  return(total_volume)
}

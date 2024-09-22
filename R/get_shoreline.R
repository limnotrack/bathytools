#' Get the shoreline from a bathymetric raster
#'
#' @inheritParams bathy_to_hypso
#' @param dTolerance numeric. The tolerance for simplifying the shoreline.
#' Default is the maximum resolution of the bathymetric raster. If set to 0, the
#' shoreline will not be simplified.
#'
#' @return sf POLYGON of the shoreline.
#' @export
#'
#' @importFrom terra clamp as.polygons
#' @importFrom sf st_as_sf st_simplify
#'
#' @examples
#' shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
#' package = "bathytools"))
#' point_data <- readRDS(system.file("extdata/depth_points.rds",
#' package = "bathytools"))
#' bathy_raster <- rasterise_bathy(shoreline = shoreline,
#'point_data = point_data, crs = 2193, res = 8)
#'shoreline2 <- get_shoreline(bathy_raster = bathy_raster)
get_shoreline <- function(bathy_raster, dTolerance = NULL) {
  if (!is(bathy_raster, "SpatRaster")) {
    stop("bathy_raster must be a SpatRaster object")
  }
  if (!is.null(dTolerance)) {
    if (!is.numeric(dTolerance)) {
      stop("dTolerance must be a numeric value")
    }
  } else {
    res <- terra::res(bathy_raster)
    dTolerance <- max(res)
  }
  r <- bathy_raster
  r <- terra::clamp(r, lower = 0, upper = Inf)
  shoreline <- r |>
    terra::as.polygons() |>
    sf::st_as_sf()
  names(shoreline) <- c("depth", "geometry")

  if (dTolerance > 0) {
    shoreline <- shoreline |>
      sf::st_simplify(dTolerance = dTolerance)
  }

  return(shoreline)
}

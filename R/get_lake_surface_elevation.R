#' Get lake surface elevation from DEM
#'
#' @inheritParams merge_bathy_dem
#'
#' @return numeric
#' @export
#'
#' @importFrom terra mask values
#'
#' @examples
#' shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
#' package = "bathytools"))
#' dem_raster <- terra::rast(system.file("extdata/dem_32m.tif",
#' package = "bathytools"))
#' lake_elev <- get_lake_surface_elevation(dem_raster = dem_raster,
#' shoreline = shoreline)

get_lake_surface_elevation <- function(dem_raster, shoreline) {
  # Check if shoreline is within the extent of dem_raster
  shoreline_in_dem <- check_contains(dem_raster, shoreline)
  if (!shoreline_in_dem) {
    stop("shoreline does not overlap with dem_raster.")
  }
  lake_dem <- terra::mask(dem_raster, shoreline)
  lake_elev <- lake_dem |>
    terra::values() |>
    median(na.rm = TRUE) |>
    round(2)
  message(paste("Lake surface elevation from DEM:", lake_elev, "m"))
  return(lake_elev)
}

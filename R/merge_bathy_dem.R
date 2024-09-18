#' Merge lake bathymetry raster with a DEM raster
#'
#' @inheritParams generate_depth_points
#' @param bathy_raster SpatRaster object with the bathymetry data.
#' @param dem_raster SpatRaster object with the DEM data.
#' @param transform_shoreline logical. If TRUE, the shoreline will be
#' transformed to match the CRS of the DEM raster.
#' @param catchment sf object of the catchment boundary. If provided, the DEM
#' raster will be cropped to the catchment boundary.
#' @param crop_dem_to_catchment logical. If TRUE, the DEM raster will be cropped
#' to the catchment boundary. Requires catchment to be provided.
#'
#' @importFrom sf st_crs st_transform st_contains st_as_sf
#' @importFrom terra resample crop mask merge
#'
#' @return SpatRaster object with the merged bathymetry and DEM data.
#' @export
#'
#' @examples
#' shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
#' package = "bathytools"))
#' catchment <- readRDS(system.file("extdata/rotoma_catchment.rds",
#'                                  package = "bathytools"))
#' point_data <- readRDS(system.file("extdata/depth_points.rds",
#' package = "bathytools"))
#' bathy_raster <- rasterise_bathy(shoreline = shoreline,
#' point_data = point_data, crs = 2193)
#' dem_raster <- terra::rast(system.file("extdata/dem_32m.tif",
#' package = "bathytools"))
#' dem_bath <- merge_bathy_dem(shoreline = shoreline, bathy_raster = bathy_raster,
#' dem_raster = dem_raster, catchment = catchment)

merge_bathy_dem <- function(shoreline, bathy_raster, dem_raster,
                            transform_shoreline = TRUE, catchment = NULL,
                            crop_dem_to_catchment = TRUE) {

  # Check if crs match
  if (sf::st_crs(dem_raster) != sf::st_crs(bathy_raster)) {
    stop("CRS of bathy_raster and dem_raster do not match.")
  }

  # Compare with CRS of the shoreline
  if (sf::st_crs(shoreline) != sf::st_crs(dem_raster)) {
    if (transform_shoreline) {
      message(strwrap("Warning: CRS of shoreline and dem_raster do
                    not match.\nReprojecting the shoreline to match the
                      dem_raster", width = 80))
      shoreline <- sf::st_transform(shoreline, sf::st_crs(dem_raster))
    }
  }

  # Check if bathy_raster is within the extent of dem_raster
  dem_sf <- spat_raster_to_sf(dem_raster)
  bathy_sf <- spat_raster_to_sf(bathy_raster)

  bathy_in_dem <- check_contains(dem_sf, bathy_sf)
  if (!bathy_in_dem) {
    stop("bathy_raster does not overlap with dem_raster.")
  }
  # Check if shoreline is within the extent of dem_raster
  shoreline_in_dem <- check_contains(dem_sf, shoreline)
  if (!shoreline_in_dem) {
    stop("shoreline does not overlap with dem_raster.")
  }
  # If catchment is provided, check if it is within the extent of dem_raster
  if (!is.null(catchment)) {
    catchment_in_dem <- check_contains(dem_sf, catchment)
    if (!catchment_in_dem) {
      stop("catchment does not overlap with dem_raster.")
    }
    if (crop_dem_to_catchment) {
      dem_raster <- dem_raster |>
        terra::crop(catchment) |>
        terra::mask(catchment)
    }
  }

  # Compare resolutions
  bathy_res <- terra::res(bathy_raster)
  dem_res <- terra::res(dem_raster)
  if (any(bathy_res != dem_res)) {
    message(strwrap(paste0("Warning: Resolutions of bathy_raster [",
                           paste0(bathy_res, collapse = ", "), "] and
                           dem_raster [",  paste0(dem_res, collapse = ", "), "]
                           do not match.\nResampling bathy_raster  to match
                           dem_raster"), width = 80))
  }
  bathy_raster <- terra::resample(bathy_raster, dem_raster)

    # Get lake elevation from DEM
  lake_elev <- get_lake_surface_elevation(dem_raster = dem_raster,
                                          shoreline = shoreline)
  # Adjust bathy to elevation
  bathy_raster <- bathy_raster + lake_elev

  # Remove lake from DEM
  dem_bath <- terra::merge(dem_raster, bathy_raster, first = FALSE)

  # Rename the variable in the raster
  names(dem_bath) <- c("elevation")

  return(dem_bath)

}

#' Check if one spatial object is contained within another
#' @param x SpatRaster or sf object
#' @param y SpatRaster or sf object
#' @return logical
#' @noRd
check_contains <- function(x, y) {
  if (is(x, "SpatRaster")) {
    x <- spat_raster_to_sf(x)
  }
  if (is(y, "SpatRaster")) {
    y <- spat_raster_to_sf(y)
  }
  all(sf::st_contains(x, y, sparse = FALSE))
}

#' Convert SpatRaster to sf polgon of
#' @param x SpatRaster
#' @return sf object
#' @noRd
spat_raster_to_sf <- function(x) {
  if (!is(x, "SpatRaster")) {
    stop("x is not a SpatRaster object")
  }
  y <- x
  y[] <- 0
  sf <- terra::as.polygons(y) |>
    sf::st_as_sf() |>
    dplyr::rename_with( .fn = ~ "raster", !dplyr::contains("geometry"))
}

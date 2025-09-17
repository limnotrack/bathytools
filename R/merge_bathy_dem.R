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
#' depth_points <- readRDS(system.file("extdata/depth_points.rds",
#' package = "bathytools"))
#' bathy_raster <- rasterise_bathy(shoreline = shoreline,
#' depth_points = depth_points, crs = 2193)
#' dem_raster <- terra::rast(system.file("extdata/dem_32m.tif",
#' package = "bathytools"))
#' dem_bath <- merge_bathy_dem(shoreline = shoreline, bathy_raster = bathy_raster,
#' dem_raster = dem_raster, catchment = catchment)

merge_bathy_dem <- function(shoreline, bathy_raster, dem_raster,
                            transform_shoreline = TRUE, catchment = NULL,
                            crop_dem_to_catchment = TRUE, blend_buffer = 10) {
  
  # ---- CRS checks ----
  dem_crs <- terra::crs(dem_raster, describe = TRUE)
  bathy_crs <- terra::crs(bathy_raster, describe = TRUE)
  if (dem_crs$code != bathy_crs$code) {
    stop("CRS of bathy_raster and dem_raster do not match.")
  }
  
  shoreline_vect <- shoreline |>
    terra::vect() |>
    terra::project(terra::crs(bathy_raster))
  
  # if (sf::st_crs(shoreline) != sf::st_crs(terra::crs(dem_raster))) {
  #   if (transform_shoreline) {
  #     message("Warning: CRS of shoreline and DEM do not match. Reprojecting shoreline.")
  #     shoreline <- sf::st_transform(shoreline, sf::st_crs(terra::crs(dem_raster)))
  #   }
  # }
  
  if (crop_dem_to_catchment) {
    if (is.null(catchment)) {
      stop("catchment must be provided when crop_dem_to_catchment = TRUE.")
    }
    dem_raster <- dem_raster |>
      terra::crop(catchment) |>
      terra::mask(catchment)
  }
  
  # ---- Resolution alignment ----
  bathy_res <- terra::res(bathy_raster)
  dem_res <- terra::res(dem_raster)
  if (!all(abs(bathy_res - dem_res) < 1e-6)) {
    message("Resolutions differ. Resampling bathy_raster to DEM resolution.")
    bathy_raster <- terra::resample(bathy_raster, dem_raster)
  }
  
  # ---- Adjust bathy to absolute elevation ----
  # rasterize lake to DEM grid (1 inside lake, NA elsewhere)
  lake_r <- terra::rasterize(shoreline_vect, dem_raster, field = 1,
                             background = NA)
  
  # distance to lake (0 on lake, positive outside)
  dist_out <- terra::distance(lake_r)
  
  # distance to land (positive inside lake): make a land mask and measure to it
  land_mask <- terra::ifel(is.na(lake_r), 1, NA)   # 1 for land, NA for lake
  dist_in_pos <- terra::distance(land_mask)        # positive inside lake (distance to land)
  dist_in <- -dist_in_pos                          # negative inside lake
  
  # combine (use dist_out where available, otherwise dist_in)
  dist_all <- terra::ifel(dist_out > 0, dist_out, dist_in)
  
  # optionally plot to check
  # terra::plot(dist_out); terra::plot(dist_in); terra::plot(dist_all)
  
  # auto blend width (e.g. twice raster resolution)
  blend_buffer <- terra::res(dist_all)[1] * 2
  
  # raw weight (SpatRaster)
  w_raw <- (dist_all + blend_buffer) / (2 * blend_buffer)
  
  # clamp to [0,1] using terra::ifel (works with SpatRaster)
  weights <- terra::ifel(w_raw < 0, 0, terra::ifel(w_raw > 1, 1, w_raw))
  
  # blend DEM and bathy_abs (bathy already converted to absolute elevation)
  # Ensure bathy_abs has values everywhere DEM has data
  bathy_abs <- bathy_raster + get_lake_surface_elevation(dem_raster, shoreline)
  bathy_full <- terra::cover(bathy_abs, dem_raster)
  
  # Now blend — no more NA outside the lake
  dem_bath <- (weights * dem_raster) + ((1 - weights) * bathy_full)
  names(dem_bath) <- "elevation"
  # terra::plot(dem_raster, main = "DEM Elevation")
  # terra::plot(dem_bath, main = "Blended Bathy and DEM")
  # terra::plot(shoreline_vect, add = TRUE)
  # terra::plet(dem_bath)
  # bathytools::plot_raster_3d(dem_bath)
  
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

#' Smooth shoreline zone of a raster
#' #' @param r SpatRaster object to smooth
#' @param shoreline sf POLYGON object of the shoreline
#' @param w numeric, size of the smoothing window. Default is 3.
#' @param fun function to apply for smoothing, default is mean.
#' @return SpatRaster object with smoothed shoreline zone
#' @noRd
smooth_shoreline <- function(r, shoreline, w = 3, fun = mean) {
  # Rasterise shoreline to match r
  shore_r <- terra::rasterize(shoreline, r, field = 1, background = NA)
  
  # Dilate shoreline mask to capture transition zone (one cell buffer)
  shore_buffer <- terra::focal(shore_r, w = matrix(1, 3, 3), fun = max, na.policy = "omit")
  shore_buffer[is.na(shore_buffer)] <- 0
  shore_buffer <- shore_buffer == 1
  
  # Apply focal smoothing to *entire raster*
  r_smooth <- terra::focal(r, w = matrix(1, w, w), fun = fun, na.policy = "omit")
  
  # Replace only shoreline zone with smoothed values
  r[shore_buffer] <- r_smooth[shore_buffer]
  
  r
}


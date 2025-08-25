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
                            crop_dem_to_catchment = TRUE) {

  # Check if crs match
  # if (terra::crs(dem_raster) != terra::crs(bathy_raster)) {
  #   stop("CRS of bathy_raster and dem_raster do not match.")
  # }

  # Compare with CRS of the shoreline
  if (sf::st_crs(shoreline) != sf::st_crs(dem_raster)) {
    if (transform_shoreline) {
      message(strwrap("Warning: CRS of shoreline and dem_raster do
                    not match.\nReprojecting the shoreline to match the
                      dem_raster", width = 80))
      shoreline <- sf::st_transform(shoreline, sf::st_crs(dem_raster))
    }
  }

  if (crop_dem_to_catchment) {
    if (is.null(catchment)) {
      stop("catchment must be provided to crop dem_raster.")
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
  if (!all(abs(bathy_res - dem_res) < 1e-6)) {
    message(strwrap(paste0("Warning: Resolutions of bathy_raster [",
                           paste0(bathy_res, collapse = ", "), "] and
                           dem_raster [",  paste0(dem_res, collapse = ", "), "]
                           do not match.\n Resampling bathy_raster  to match
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
  
  # Smooth shoreline
  dem_bath_smooth <- smooth_shoreline(dem_bath, shoreline)
  
  mm <- terra::minmax(bathy_raster)
  lake_depth <- round(mm[2] - mm[1], 2)
  ext_hyps <- dem_to_hypsograph(shoreline = shoreline, dem_bath = dem_bath, 
                                lake_depth = lake_depth,
                                ext_elev = 3)
  
  # Rename the variable in the raster
  names(dem_bath) <- c("elevation")

  return(dem_bath)

}

#' Merge lake bathymetry raster with a DEM raster (smooth blend at shoreline)
#'
#' @inheritParams generate_depth_points
#' @param bathy_raster SpatRaster object with bathymetry depths (negative values).
#' @param dem_raster SpatRaster object with DEM elevations.
#' @param transform_shoreline logical. If TRUE, reprojects shoreline to match DEM CRS.
#' @param catchment sf object of the catchment boundary.
#' @param crop_dem_to_catchment logical. If TRUE, crops DEM to catchment boundary.
#' @param blend_buffer numeric. Width (map units) of the shoreline blending zone.
#'
#' @importFrom sf st_crs st_transform st_as_sf st_bbox st_sfc
#' @importFrom terra crs crop mask resample merge cover distance ifel ext
#'
#' @return SpatRaster object with merged bathymetry and DEM data.
#' @export
merge_bathy_dem <- function(shoreline, bathy_raster, dem_raster,
                            transform_shoreline = TRUE, catchment = NULL,
                            crop_dem_to_catchment = TRUE, blend_buffer = 10) {
  
  # ---- CRS checks ----
  if (terra::crs(dem_raster) != terra::crs(bathy_raster)) {
    stop("CRS of bathy_raster and dem_raster do not match.")
  }
  
  if (sf::st_crs(shoreline) != sf::st_crs(terra::crs(dem_raster))) {
    if (transform_shoreline) {
      message("Warning: CRS of shoreline and DEM do not match. Reprojecting shoreline.")
      shoreline <- sf::st_transform(shoreline, sf::st_crs(terra::crs(dem_raster)))
    }
  }
  
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
  lake_r <- terra::rasterize(terra::vect(shoreline), dem_raster, field = 1,
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
  terra::plot(dem_raster, main = "DEM Elevation")
  terra::plot(dem_bath, main = "Blended Bathy and DEM")
  
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


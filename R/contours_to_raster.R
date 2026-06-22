#' Convert lake depth contours to a raster
#'
#' @param contours An sf object containing depth contour lines or polygons,
#'   with a numeric depth attribute.
#' @param depth_col Character. Name of the column in `contours` containing
#'   depth values.
#' @param res Numeric. Output raster resolution in the CRS units (e.g. metres).
#' @param crs An sf/terra-compatible CRS string, EPSG integer, or NULL to
#'   inherit from `contours`.
#' @param method Character. Interpolation method passed to
#'   \code{terra::interpIDW} or \code{"near"} for nearest-contour assignment.
#'   Use \code{"idw"} (default) for smooth interpolation.
#' @param idw_power Numeric. IDW power parameter (default 2).
#'
#' @return A \code{terra::SpatRaster} of interpolated depth values, masked to
#'   the lake boundary derived from the outermost contour.
#'   
#'
#' @examples
#' \dontrun{
#'   contours <- sf::st_read("lake_contours.gpkg")
#'   depth_raster <- contours_to_raster(contours, depth_col = "depth_m", res = 5)
#'   terra::plot(depth_raster)
#' }
contours_to_raster <- function(shoreline,
                               contours,
                               depth_col  = "depth",
                               res        = 10,
                               crs        = NULL,
                               method     = c("idw", "near"),
                               idw_power  = 2) {
  
  method <- match.arg(method)
  
  # ── 1. Validate inputs ───────────────────────────────────────────────────────
  if (!inherits(contours, "sf")) {
    stop("`contours` must be an sf object.")
  }
  if (!depth_col %in% names(contours)) {
    stop(sprintf("Column '%s' not found in `contours`.", depth_col))
  }
  if (!is.numeric(contours[[depth_col]])) {
    stop(sprintf("Column '%s' must be numeric.", depth_col))
  }
  
  # ── 2. Re-project if requested ───────────────────────────────────────────────
  if (!is.null(crs)) {
    contours <- sf::st_transform(contours, crs)
  }
  
  # Add shoreline to contours
  shoreline_vect <- terra::vect(shoreline)
  shoreline[[depth_col]] <- 0   # add depth=0 to shoreline for interpolation  
  contours <- dplyr::bind_rows(shoreline, contours)
  
  # ── 3. Ensure geometries are LINESTRING / MULTILINESTRING ────────────────────
  #       (convert polygon rings to lines so rasterisation works uniformly)
  geom_type <- unique(sf::st_geometry_type(contours, by_geometry = FALSE))
  
  if (any(grepl(c("POLYGON|GEOMETRY"), geom_type, ignore.case = TRUE))) {
    contours <- sf::st_cast(contours, "MULTILINESTRING")
  }
  
  # ── 4. Build lake boundary mask from the outermost (shallowest) contour ──────
  #       Union all lines → convex-hull polygon used as the lake extent mask.
  lake_boundary <- contours |>
    sf::st_union() |>
    sf::st_convex_hull() |>
    sf::st_as_sf()
  
  # ── 5. Create template raster aligned to the contour bbox ───────────────────
  bbox   <- terra::ext(terra::vect(contours))
  template_rast <- terra::rast(
    extent = bbox,
    resolution = res,
    crs = terra::crs(terra::vect(contours))
  )
  
  # ── 6. Rasterise contour lines (burn depth values) ──────────────────────────
  contours_vect <- terra::vect(contours)
  
  contour_rast <- terra::rasterize(
    x     = contours_vect,
    y     = template_rast,
    field = depth_col,
    fun   = "mean"      # average where contours overlap a cell
  )
  
  # ── 7. Interpolate across the full lake extent ───────────────────────────────
  depth_rast <- if (method == "idw") {
    
    # Convert rasterised contour cells to points for IDW
    contour_pts <- terra::as.points(contour_rast, na.rm = TRUE)
    
    # Estimate max gap by sampling a subset of contour point coordinates
    # and finding the maximum nearest-neighbour distance among them.
    coords     <- terra::crds(contour_pts)
    n_sample   <- min(500L, nrow(coords))
    idx        <- sample(nrow(coords), n_sample)
    sample_pts <- coords[idx, , drop = FALSE]
    
    # For each sampled point find its nearest neighbour distance
    nn_dists <- apply(sample_pts, 1, \(p) {
      dists <- sqrt((sample_pts[, 1] - p[1])^2 + (sample_pts[, 2] - p[2])^2)
      sort(dists)[2]   # [1] is self (0), [2] is nearest neighbour
    })
    
    safe_radius <- max(nn_dists) * 2   # double the max gap for safety
    
    terra::interpIDW(
      x         = template_rast,
      y         = contour_pts,
      field     = depth_col,
      radius    = safe_radius,
      power     = idw_power,
      maxPoints = 12
    )
    
  } else {
    
    # Nearest-contour assignment via focal distance fill
    terra::distance(contour_rast) |>   # distance raster (helper — not used directly)
      {\(.) terra::focal(
        contour_rast,
        w   = terra::focalMat(contour_rast, res * 5, type = "circle"),
        fun = "mean",
        na.policy = "only",
        na.rm = TRUE
      )}()
  }
  
  # ── 8. Mask to lake boundary ─────────────────────────────────────────────────
  lake_vect  <- terra::vect(shoreline)
  depth_rast <- terra::mask(depth_rast, lake_vect)
  terra::plot(depth_rast)
  
  
  names(depth_rast) <- depth_col
  depth_rast
}

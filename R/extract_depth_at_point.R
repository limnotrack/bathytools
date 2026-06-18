#' Extract Lake Depth at a Given Point
#'
#' This function extracts the lake depth at a specific point location from
#' various data sources: bathymetry raster, point depth data, or contour data.
#' If multiple data sources are provided, the function will use them in the
#' following priority order: bathymetry raster > point data > contours.
#'
#' @param x numeric or sf POINT. Either a numeric vector of length 2 containing
#'   the x and y coordinates (c(x, y)), or an sf POINT object. If numeric
#'   coordinates are provided, the `crs` parameter must also be specified.
#' @param bathy_raster SpatRaster object with the bathymetry data. Optional if
#'   `depth_points` or `contours` is provided. If provided along with other
#'   data sources, this will be used preferentially.
#' @param depth_points sf POINT object with depth data. Must contain a 'depth'
#'   column. Optional if `bathy_raster` or `contours` is provided. If provided
#'   along with `contours` (but not `bathy_raster`), this will be used preferentially.
#' @param contours sf LINESTRING or MULTILINESTRING object with contour data.
#'   Must contain a 'depth' column. Optional if `bathy_raster` or `depth_points`
#'   is provided. This has the lowest priority if multiple data sources are provided.
#' @param crs numeric or character. Coordinate reference system (CRS) of the
#'   input coordinates if `x` is a numeric vector. Can be an EPSG code or proj4
#'   string. Required if `x` is numeric; ignored if `x` is an sf object.
#' @param method character. Method to use for extraction. Options are:
#'   \itemize{
#'     \item "bilinear": bilinear interpolation from raster (default for raster)
#'     \item "simple": nearest cell value from raster
#'     \item "nearest": nearest neighbor from point data or contours (default for point/contour data)
#'     \item "idw": inverse distance weighting from point data
#'   }
#'   Default is "bilinear" for raster data and "nearest" for point/contour data.
#' @param max_dist numeric. Maximum distance (in map units) to search for
#'   nearest neighbor or for IDW interpolation when using point or contour data.
#'   Default is Inf (no distance limit).
#' @param idw_power numeric. Power parameter for inverse distance weighting.
#'   Default is 2.
#' @param n_neighbors integer. Number of nearest neighbors to use for IDW
#'   interpolation. Default is 5.
#'
#' @return numeric. The depth value at the specified point. Returns NA if no
#'   valid depth can be extracted (e.g., point is outside the raster extent,
#'   or no nearby points/contours are found within max_dist).
#'
#' @export
#'
#' @importFrom terra extract crs vect
#' @importFrom sf st_as_sf st_coordinates st_crs st_transform st_distance st_nearest_feature st_geometry_type
#' @importFrom units drop_units
#'
#' @examples
#' \dontrun{
#' # Load example data
#' shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
#'                                  package = "bathytools"))
#' depth_points <- readRDS(system.file("extdata/depth_points.rds",
#'                                   package = "bathytools"))
#'
#' # Create bathymetry raster
#' bathy <- rasterise_bathy(shoreline = shoreline, depth_points = depth_points,
#'                          crs = 2193, res = 8)
#'
#' # Extract depth at a point from raster
#' depth <- extract_depth_at_point(x = c(2823700, 6404300),
#'                                 bathy_raster = bathy,
#'                                 crs = 2193)
#'
#' # Extract depth from point data
#' depth <- extract_depth_at_point(x = c(2823700, 6404300),
#'                                 depth_points = depth_points,
#'                                 crs = 2193,
#'                                 method = "nearest")
#'
#' # Extract depth using IDW from point data
#' depth <- extract_depth_at_point(x = c(2823700, 6404300),
#'                                 depth_points = depth_points,
#'                                 crs = 2193,
#'                                 method = "idw",
#'                                 n_neighbors = 5)
#'
#' # Extract depth from contours
#' contours <- get_contours(bathy_raster = bathy)
#' depth <- extract_depth_at_point(x = c(2823700, 6404300),
#'                                 contours = contours,
#'                                 crs = 2193,
#'                                 method = "nearest")
#' }
#'
extract_depth_at_point <- function(x,
                                   bathy_raster = NULL,
                                   depth_points = NULL,
                                   contours = NULL,
                                   crs = NULL,
                                   method = NULL,
                                   max_dist = Inf,
                                   idw_power = 2,
                                   n_neighbors = 5) {
  
  # Input validation
  if (is.null(bathy_raster) && is.null(depth_points) && is.null(contours)) {
    stop("At least one of 'bathy_raster', 'depth_points', or 'contours' must be provided.")
  }
  
  # Convert x to sf point if it's numeric
  if (is.numeric(x)) {
    if (length(x) != 2) {
      stop("If 'x' is numeric, it must be a vector of length 2 (x, y coordinates).")
    }
    if (is.null(crs)) {
      stop("If 'x' is numeric, 'crs' must be specified.")
    }
    # Create sf point from coordinates
    point_sf <- sf::st_as_sf(
      data.frame(x = x[1], y = x[2]),
      coords = c("x", "y"),
      crs = crs
    )
  } else if (inherits(x, "sf") || inherits(x, "sfc")) {
    point_sf <- x
    # Check that it's a POINT geometry
    geom_type <- unique(as.character(sf::st_geometry_type(point_sf)))
    if (length(geom_type) != 1 || geom_type != "POINT") {
      stop("If 'x' is an sf object, it must be a POINT geometry.")
    }
  } else {
    stop("'x' must be either a numeric vector of length 2 or an sf POINT object.")
  }
  
  # Ensure point is a single feature
  if (nrow(point_sf) > 1) {
    warning("Multiple points provided. Only the first point will be used.")
    point_sf <- point_sf[1, ]
  }
  
  # Determine which data source to use (priority: raster > points > contours)
  if (!is.null(bathy_raster)) {
    # Extract from raster
    if (!inherits(bathy_raster, "SpatRaster")) {
      stop("'bathy_raster' must be a SpatRaster object.")
    }
    
    # Set default method for raster
    if (is.null(method)) {
      method <- "bilinear"
    }
    
    # Validate method for raster
    if (!method %in% c("bilinear", "simple")) {
      stop("For raster extraction, method must be 'bilinear' or 'simple'.")
    }
    
    # Transform point to raster CRS
    raster_crs <- terra::crs(bathy_raster)
    point_transformed <- sf::st_transform(point_sf, sf::st_crs(raster_crs))
    
    # Convert to terra vect for extraction
    point_vect <- terra::vect(point_transformed)
    
    # Extract value
    depth_value <- terra::extract(bathy_raster, point_vect, method = method)
    
    # Return the depth value (extract returns a data.frame)
    return(as.numeric(depth_value[1, 2]))
    
  } else if (!is.null(depth_points)) {
    # Extract from point data
    if (!inherits(depth_points, "sf")) {
      stop("'depth_points' must be an sf object.")
    }
    
    if (!"depth" %in% names(depth_points)) {
      stop("'depth_points' must contain a 'depth' column.")
    }
    
    # Set default method for point data
    if (is.null(method)) {
      method <- "nearest"
    }
    
    # Validate method for point data
    if (!method %in% c("nearest", "idw")) {
      stop("For point data extraction, method must be 'nearest' or 'idw'.")
    }
    
    # Transform point to depth_points CRS
    point_transformed <- sf::st_transform(point_sf, sf::st_crs(depth_points))
    
    if (method == "nearest") {
      # Find nearest point
      nearest_idx <- sf::st_nearest_feature(point_transformed, depth_points)
      
      # Calculate distance to nearest point
      dist_to_nearest <- sf::st_distance(point_transformed, depth_points[nearest_idx, ])
      dist_to_nearest <- units::drop_units(dist_to_nearest)[1, 1]
      
      # Check if within max_dist
      if (dist_to_nearest > max_dist) {
        return(NA_real_)
      }
      
      # Return depth of nearest point
      return(depth_points$depth[nearest_idx])
      
    } else if (method == "idw") {
      # Calculate distances to all points
      distances <- sf::st_distance(point_transformed, depth_points)
      distances <- units::drop_units(distances)[1, ]
      
      # Filter by max_dist
      within_dist <- distances <= max_dist
      if (sum(within_dist) == 0) {
        return(NA_real_)
      }
      
      # Get n nearest neighbors
      n_use <- min(n_neighbors, sum(within_dist))
      nearest_indices <- order(distances)[1:n_use]
      nearest_distances <- distances[nearest_indices]
      
      # Handle case where point coincides with a data point
      if (any(nearest_distances == 0)) {
        zero_idx <- which(nearest_distances == 0)[1]
        return(depth_points$depth[nearest_indices[zero_idx]])
      }
      
      # Calculate IDW weights
      weights <- 1 / (nearest_distances ^ idw_power)
      weights <- weights / sum(weights)
      
      # Calculate weighted depth
      depth_value <- sum(depth_points$depth[nearest_indices] * weights)
      return(depth_value)
    }
    
  } else if (!is.null(contours)) {
    # Extract from contour data
    if (!inherits(contours, "sf")) {
      stop("'contours' must be an sf object.")
    }
    
    if (!"depth" %in% names(contours)) {
      stop("'contours' must contain a 'depth' column.")
    }
    
    # Set default method for contour data
    if (is.null(method)) {
      method <- "nearest"
    }
    
    # Validate method for contour data
    if (method != "nearest") {
      stop("For contour data extraction, only method 'nearest' is supported.")
    }
    
    # Transform point to contours CRS
    point_transformed <- sf::st_transform(point_sf, sf::st_crs(contours))
    
    # Find nearest contour
    nearest_idx <- sf::st_nearest_feature(point_transformed, contours)
    
    # Calculate distance to nearest contour
    dist_to_nearest <- sf::st_distance(point_transformed, contours[nearest_idx, ])
    dist_to_nearest <- units::drop_units(dist_to_nearest)[1, 1]
    
    # Check if within max_dist
    if (dist_to_nearest > max_dist) {
      return(NA_real_)
    }
    
    # Return depth of nearest contour
    return(contours$depth[nearest_idx])
  }
}

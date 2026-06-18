#' Extract Lake Depth at a Given Point
#'
#' This function extracts the lake depth at a specific point location from
#' various data sources: bathymetry raster, point depth data, or contour data.
#' If multiple data sources are provided, the function will use them in the
#' following priority order: bathymetry raster > point data > contours.
#'
#' The function performs boundary checks to verify that the point is within the
#' extent of the data source before attempting extraction. For raster data, the
#' function will return NA with a warning if the point is outside the raster extent.
#' For point and contour data, the function will issue a warning but continue with
#' the search, as points may still be found within the specified max_dist.
#'
#' The function uses the cli package to provide informative messages about the
#' extraction method, boundary checks, and results, making it easier to understand
#' what the function is doing.
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
#' @importFrom terra extract crs vect ext
#' @importFrom sf st_as_sf st_coordinates st_crs st_transform st_distance st_nearest_feature st_geometry_type st_bbox
#' @importFrom units drop_units
#' @importFrom cli cli_inform cli_warn cli_abort
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
    cli::cli_abort("At least one of 'bathy_raster', 'depth_points', or 'contours' must be provided.")
  }
  
  # Convert x to sf point if it's numeric
  if (is.numeric(x)) {
    if (length(x) != 2) {
      cli::cli_abort("If 'x' is numeric, it must be a vector of length 2 (x, y coordinates).")
    }
    if (is.null(crs)) {
      cli::cli_abort("If 'x' is numeric, 'crs' must be specified.")
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
      cli::cli_abort("If 'x' is an sf object, it must be a POINT geometry.")
    }
  } else {
    cli::cli_abort("'x' must be either a numeric vector of length 2 or an sf POINT object.")
  }
  
  # Ensure point is a single feature
  if (nrow(point_sf) > 1) {
    cli::cli_warn(c("!" = "Multiple points provided. Only the first point will be used."))
    point_sf <- point_sf[1, ]
  }
  
  # Determine which data source to use (priority: raster > points > contours)
  if (!is.null(bathy_raster)) {
    # Extract from raster
    if (!inherits(bathy_raster, "SpatRaster")) {
      cli::cli_abort("'bathy_raster' must be a SpatRaster object.")
    }
    
    # Set default method for raster
    if (is.null(method)) {
      method <- "bilinear"
    }
    
    # Validate method for raster
    if (!method %in% c("bilinear", "simple")) {
      cli::cli_abort("For raster extraction, method must be 'bilinear' or 'simple'.")
    }
    
    cli::cli_inform(c("i" = "Extracting depth from bathymetry raster using {.field {method}} method."))
    
    # Transform point to raster CRS
    raster_crs <- terra::crs(bathy_raster)
    point_transformed <- sf::st_transform(point_sf, sf::st_crs(raster_crs))
    
    # Check if point is within raster extent
    raster_ext <- terra::ext(bathy_raster)
    point_coords <- sf::st_coordinates(point_transformed)
    
    if (point_coords[1] < raster_ext[1] || point_coords[1] > raster_ext[2] ||
        point_coords[2] < raster_ext[3] || point_coords[2] > raster_ext[4]) {
      cli::cli_warn(c("!" = "Point is outside the bathymetry raster extent.",
                      "i" = "Returning NA."))
      return(NA_real_)
    }
    
    cli::cli_inform(c("v" = "Point is within raster extent."))
    
    # Convert to terra vect for extraction
    point_vect <- terra::vect(point_transformed)
    
    # Extract value
    depth_value <- terra::extract(bathy_raster, point_vect, method = method)
    
    # Return the depth value (extract returns a data.frame)
    return(as.numeric(depth_value[1, 2]))
    
  } else if (!is.null(depth_points)) {
    # Extract from point data
    if (!inherits(depth_points, "sf")) {
      cli::cli_abort("'depth_points' must be an sf object.")
    }
    
    if (!"depth" %in% names(depth_points)) {
      cli::cli_abort("'depth_points' must contain a 'depth' column.")
    }
    
    # Set default method for point data
    if (is.null(method)) {
      method <- "nearest"
    }
    
    # Validate method for point data
    if (!method %in% c("nearest", "idw")) {
      cli::cli_abort("For point data extraction, method must be 'nearest' or 'idw'.")
    }
    
    cli::cli_inform(c("i" = "Extracting depth from point data using {.field {method}} method."))
    
    # Transform point to depth_points CRS
    point_transformed <- sf::st_transform(point_sf, sf::st_crs(depth_points))
    
    # Check if point is within extent of depth points
    depth_bbox <- sf::st_bbox(depth_points)
    point_coords <- sf::st_coordinates(point_transformed)
    
    if (point_coords[1] < depth_bbox["xmin"] || point_coords[1] > depth_bbox["xmax"] ||
        point_coords[2] < depth_bbox["ymin"] || point_coords[2] > depth_bbox["ymax"]) {
      cli::cli_warn(c("!" = "Point is outside the extent of depth points.",
                      "i" = "Continuing with {method} search (may return NA if beyond max_dist)."))
    } else {
      cli::cli_inform(c("v" = "Point is within depth points extent."))
    }
    
    if (method == "nearest") {
      # Find nearest point
      nearest_idx <- sf::st_nearest_feature(point_transformed, depth_points)
      
      # Calculate distance to nearest point
      dist_to_nearest <- sf::st_distance(point_transformed, depth_points[nearest_idx, ])
      dist_to_nearest <- units::drop_units(dist_to_nearest)[1, 1]
      
      # Check if within max_dist
      if (dist_to_nearest > max_dist) {
        cli::cli_warn(c("!" = "Nearest point is {round(dist_to_nearest, 2)} units away, exceeding max_dist ({max_dist}).",
                        "i" = "Returning NA."))
        return(NA_real_)
      }
      
      cli::cli_inform(c("v" = "Found nearest point at distance {round(dist_to_nearest, 2)} units."))
      
      # Return depth of nearest point
      return(depth_points$depth[nearest_idx])
      
    } else if (method == "idw") {
      # Calculate distances to all points
      distances <- sf::st_distance(point_transformed, depth_points)
      distances <- units::drop_units(distances)[1, ]
      
      # Filter by max_dist
      within_dist <- distances <= max_dist
      if (sum(within_dist) == 0) {
        cli::cli_warn(c("!" = "No depth points found within max_dist ({max_dist}).",
                        "i" = "Returning NA."))
        return(NA_real_)
      }
      
      # Get n nearest neighbors
      n_use <- min(n_neighbors, sum(within_dist))
      
      # Ensure we have at least one neighbor
      if (n_use == 0) {
        cli::cli_warn(c("!" = "No neighbors available for IDW interpolation.",
                        "i" = "Returning NA."))
        return(NA_real_)
      }
      
      cli::cli_inform(c("i" = "Using {n_use} nearest neighbor{?s} for IDW interpolation."))
      
      nearest_indices <- order(distances)[1:n_use]
      nearest_distances <- distances[nearest_indices]
      
      # Handle case where point coincides with a data point
      if (any(nearest_distances == 0)) {
        zero_idx <- which(nearest_distances == 0)[1]
        cli::cli_inform(c("v" = "Point coincides exactly with a depth point.",
                          "i" = "Returning exact depth value."))
        return(depth_points$depth[nearest_indices[zero_idx]])
      }
      
      # Calculate IDW weights
      weights <- 1 / (nearest_distances ^ idw_power)
      weights <- weights / sum(weights)
      
      # Calculate weighted depth
      depth_value <- sum(depth_points$depth[nearest_indices] * weights)
      cli::cli_inform(c("v" = "IDW interpolation complete."))
      return(depth_value)
    }
    
  } else if (!is.null(contours)) {
    # Extract from contour data
    if (!inherits(contours, "sf")) {
      cli::cli_abort("'contours' must be an sf object.")
    }
    
    if (!"depth" %in% names(contours)) {
      cli::cli_abort("'contours' must contain a 'depth' column.")
    }
    
    # Set default method for contour data
    if (is.null(method)) {
      method <- "nearest"
    }
    
    # Validate method for contour data
    if (method != "nearest") {
      cli::cli_abort("For contour data extraction, only method 'nearest' is supported.")
    }
    
    cli::cli_inform(c("i" = "Extracting depth from contour data using {.field {method}} method."))
    
    # Transform point to contours CRS
    point_transformed <- sf::st_transform(point_sf, sf::st_crs(contours))
    
    # Check if point is within extent of contours
    contours_bbox <- sf::st_bbox(contours)
    point_coords <- sf::st_coordinates(point_transformed)
    
    if (point_coords[1] < contours_bbox["xmin"] || point_coords[1] > contours_bbox["xmax"] ||
        point_coords[2] < contours_bbox["ymin"] || point_coords[2] > contours_bbox["ymax"]) {
      cli::cli_warn(c("!" = "Point is outside the extent of contours.",
                      "i" = "Continuing with {method} search (may return NA if beyond max_dist)."))
    } else {
      cli::cli_inform(c("v" = "Point is within contours extent."))
    }
    
    # Find nearest contour
    nearest_idx <- sf::st_nearest_feature(point_transformed, contours)
    
    # Calculate distance to nearest contour
    dist_to_nearest <- sf::st_distance(point_transformed, contours[nearest_idx, ])
    dist_to_nearest <- units::drop_units(dist_to_nearest)[1, 1]
    
    # Check if within max_dist
    if (dist_to_nearest > max_dist) {
      cli::cli_warn(c("!" = "Nearest contour is {round(dist_to_nearest, 2)} units away, exceeding max_dist ({max_dist}).",
                      "i" = "Returning NA."))
      return(NA_real_)
    }
    
    cli::cli_inform(c("v" = "Found nearest contour at distance {round(dist_to_nearest, 2)} units."))
    
    # Return depth of nearest contour
    return(contours$depth[nearest_idx])
  }
}

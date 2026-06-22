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
#' For point and contour data, the function will issue a warning if the point is
#' outside the data extent but will continue searching, as features may still be
#' found within the specified max_dist (though NA will be returned if none are found).
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
#' @param shoreline sf POLYGON or MULTIPOLYGON object representing the 
#' shoreline. Optional but required if `contours` are provided, as the shoreline 
#' is used to define the lake boundary for interpolation. Ignored if `contours` 
#' is not provided.
#' @param islands sf POLYGON or MULTIPOLYGON object representing any islands in
#'  the lake. Default is NULL. Optional but recommended if islands are present, 
#'  as they can affect the interpolation of depth points and contours.
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
#'     \item "nearest": nearest neighbor from point data or contours
#'     \item "idw": inverse distance weighting from point data
#'   }
#'   Default is "bilinear" for raster data and "nearest" for point/contour data.
#' @param max_dist numeric. Maximum distance (in map units) to search for
#'   nearest neighbor or for IDW interpolation when using point or contour data.
#'   Default is Inf (no distance limit).
#' @param idw_power numeric. Power parameter for inverse distance weighting.
#'   Default is 2.
#' @param n_neighbors integer. Number of nearest neighbors to use for IDW
#'   interpolation. Default is 50.
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
#' @importFrom rlang arg_match
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
                                   shoreline = NULL,
                                   islands = NULL,
                                   contours = NULL,
                                   crs = NULL,
                                   method = c("bilinear", "simple", "nearest", "idw"),
                                   max_dist = Inf,
                                   idw_power = 2,
                                   n_neighbors = 50) {
  
  # Input validation
  if (is.null(bathy_raster) && is.null(depth_points) && is.null(contours)) {
    cli::cli_abort("At least one of 'bathy_raster', 'depth_points', or 'contours' must be provided.")
  }
  method <- rlang::arg_match(method)
  
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
  
  # If contours supplied, convert to depth points
  if (!is.null(shoreline) && !is.null(contours) && method != "nearest") {
    missing_islands <- detect_islands(shoreline)
    if (!is.null(missing_islands)) {
      if (is.null(islands)) {
        islands <- missing_islands
      } else {
        islands <- dplyr::bind_rows(islands, missing_islands)
      }
    }
    shoreline <- detect_shoreline(shoreline)
    res <- est_res(shoreline)
    depth_points <- generate_depth_points(shoreline = shoreline, 
                                          islands = islands,
                                          contours = contours, 
                                          res = res)
    
    bathy_raster <- interpolate_points(depth_points = depth_points, 
                                       shoreline = shoreline,
                                       islands = islands,
                                       res = res,
                                       print_plot = FALSE)
    
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
    
  } else if (!is.null(contours) || !is.null(depth_points)) {
    
    # --- Normalise whichever source was supplied ----------------------------
    
    if (!is.null(contours)) {
      
      if (!inherits(contours, "sf")) {
        cli::cli_abort("'contours' must be an sf object.")
      }
      
      if (!"depth" %in% names(contours)) {
        cli::cli_abort("'contours' must contain a 'depth' column.")
      }
      
      source_data   <- contours
      source_label  <- "contours"
      
    } else {
      
      # Convert plain data frame to sf
      if (is.data.frame(depth_points) && !inherits(depth_points, "sf")) {
        lon_lat_cols <- guess_lonlat_cols(depth_points)
        if (is.null(crs)) crs <- 4326
        depth_points <- sf::st_as_sf(depth_points, coords = lon_lat_cols, crs = crs)
      }
      
      if (!inherits(depth_points, "sf")) {
        cli::cli_abort("'depth_points' must be an sf object.")
      }
      
      if (!"depth" %in% names(depth_points)) {
        cli::cli_abort("'depth_points' must contain a 'depth' column.")
      }
      
      source_data   <- depth_points
      source_label  <- "depth_points"
      
    }
    
    # --- Method validation --------------------------------------------------
    
    if (is.null(method)) {
      method <- "idw"
    }
    
    if (!method %in% c("nearest", "idw")) {
      cli::cli_abort("For point/contour data extraction, method must be 'nearest' or 'idw'.")
    }
    
    cli::cli_inform(c("i" = "Extracting depth from {.field {source_label}} using {.field {method}} method."))
    
    # --- Transform & extent check -------------------------------------------
    
    point_transformed <- sf::st_transform(point_sf, sf::st_crs(source_data))
    
    depth_bbox   <- sf::st_bbox(source_data)
    point_coords <- sf::st_coordinates(point_transformed)
    
    if (point_coords[1] < depth_bbox["xmin"] || point_coords[1] > depth_bbox["xmax"] ||
        point_coords[2] < depth_bbox["ymin"] || point_coords[2] > depth_bbox["ymax"]) {
      cli::cli_warn(c(
        "!" = "Point is outside the extent of {.field {source_label}}.",
        "i" = "Continuing with {.field {method}} search (may return NA if beyond max_dist)."
      ))
    } else {
      cli::cli_inform(c("v" = "Point is within {.field {source_label}} extent."))
    }
    
    # --- Nearest ------------------------------------------------------------
    
    if (method == "nearest") {
      
      nearest_idx     <- sf::st_nearest_feature(point_transformed, source_data)
      dist_to_nearest <- sf::st_distance(point_transformed, source_data[nearest_idx, ])
      dist_to_nearest <- units::drop_units(dist_to_nearest)[1, 1]
      
      if (dist_to_nearest > max_dist) {
        cli::cli_warn(c(
          "!" = "Nearest feature is {round(dist_to_nearest, 2)} units away, exceeding max_dist ({max_dist}).",
          "i" = "Returning NA."
        ))
        return(NA_real_)
      }
      
      cli::cli_inform(c("v" = "Found nearest feature at distance {round(dist_to_nearest, 2)} units."))
      return(source_data$depth[nearest_idx])
      
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
      n_use <- pmin(n_neighbors, sum(within_dist))
      
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
      
      tm_shape(depth_points[nearest_indices, ]) +
        tm_dots(col = "depth", size = 0.5) +
        tm_shape(point_transformed) +
        tm_dots(col = "red", size = 0.5) +
        tm_layout(main.title = "IDW Interpolation Neighbors")
      
      # Calculate weighted depth
      depth_value <- sum(depth_points$depth[nearest_indices] * weights)
      cli::cli_inform(c("v" = "IDW interpolation complete."))
      return(depth_value)
    }
  }
}

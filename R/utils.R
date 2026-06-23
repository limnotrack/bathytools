#' Guess longitude and latitude columns in a data frame.
#' This function looks for common column names that might represent longitude and latitude.
#' It returns a named vector with the guessed longitude and latitude column names.
#' If it cannot find suitable columns, it throws an error with the available column names.
#' @param df A data frame to inspect for longitude and latitude columns.
#' @return A named vector with 'lon' and 'lat' as names and the
#' guessed column names as values.
#' @noRd
#' @importFrom cli cli_abort
guess_lonlat_cols <- function(df) {
  nms <- names(df)
  lon_candidates <- c("lon", "lng", "long", "longitude", "x")
  lat_candidates <- c("lat", "latitude", "y")
  
  lon_col <- nms[match(TRUE, tolower(nms) %in% lon_candidates)]
  lat_col <- nms[match(TRUE, tolower(nms) %in% lat_candidates)]
  
  if (is.na(lon_col) || is.na(lat_col)) {
    cli::cli_abort(c(
      "x" = "Could not guess longitude and latitude columns.",
      "i" = "Found columns: {paste(nms, collapse = ', ')}"
    ))
  }
  
  c(lon = lon_col, lat = lat_col)
}

#' Estimate an appropriate raster resolution for a spatial polygon
#'
#' Scales resolution based on polygon area using a power-law relationship,
#' so small areas get fine resolution and large areas get coarser resolution.
#'
#' @param x   An sf polygon or multipolygon object
#' @param target_cells  Target number of cells across the shorter dimension (default 500)
#' @param min_res  Minimum resolution in CRS units (default NULL = no floor)
#' @param max_res  Maximum resolution in CRS units (default NULL = no ceiling)
#' @param snap_to  Optional numeric vector of "nice" values to snap to (e.g. c(5, 10, 25, 50, 100))
#'
#' @return A single numeric value: the estimated resolution in CRS units
#'
#' @examples
#' est_res(nz_coast)                        # auto-scale
#' est_res(nz_coast, target_cells = 200)    # coarser
#' est_res(nz_coast, snap_to = c(10, 25, 50, 100, 250, 500))

est_res <- function(x,
                    target_cells = 500,
                    min_res      = NULL,
                    max_res      = NULL,
                    snap_to      = c(2, 4, 8, 16, 20, 24, 32, 40, 64)) {
  
  stopifnot(inherits(x, "sf") || inherits(x, "sfc"))
  
  # Warn if geometry is not in metric CRS (resolution will be in CRS units, which may not be metres)
  if (!grepl("metre|meter|foot|us-ft|crs\\ =\\s*\\d+", sf::st_crs(x)$units_gdal, 
             ignore.case = TRUE)) {
    cli::cli_warn(c(
      "x" = "Geometry does not appear to be in a metric CRS. Estimated resolution will be in CRS units, which may not be metres.",
      "i" = "CRS units: {sf::st_crs(x)$units_gdal}"
    ))
  }
  
  bb     <- sf::st_bbox(x)
  width  <- as.numeric(bb["xmax"] - bb["xmin"])
  height <- as.numeric(bb["ymax"] - bb["ymin"])
  
  # Drive resolution off the shorter side so cells stay roughly square
  shorter <- min(width, height)
  res     <- shorter / target_cells
  
  # Apply floor / ceiling if requested
  if (!is.null(min_res)) res <- max(res, min_res)
  if (!is.null(max_res)) res <- min(res, max_res)
  
  # Snap to nearest "nice" value if a vector is provided
  if (!is.null(snap_to)) {
    diffs <- abs(snap_to - res)
    res   <- snap_to[which.min(diffs)]
  }
  
  return(res)
}

#' Detect islands within a lake multipolygon
#'
#' Islands are encoded as inner rings (holes) in the exterior polygon.
#' This function extracts those holes and returns them as proper polygons.
#'
#' @param x An sf object with POLYGON or MULTIPOLYGON geometry
#' @return An sf object containing island polygons, or NULL if none found
#' @export
#' @importFrom sf st_geometry st_polygon st_multipolygon st_sfc st_crs st_sf
detect_islands <- function(x) {
  stopifnot(inherits(x, "sf") || inherits(x, "sfc"))
  
  geom <- if (inherits(x, "sf")) sf::st_geometry(x) else x
  
  # Collect all inner rings across every (multi)polygon
  island_rings <- list()
  
  for (i in seq_along(geom)) {
    g <- geom[[i]]
    
    # Normalise: treat both POLYGON and MULTIPOLYGON uniformly
    polys <- if (inherits(g, "MULTIPOLYGON")) {
      lapply(seq_along(g), function(j) g[[j]])
    } else {
      list(g)
    }
    
    for (poly in polys) {
      # poly is a list of rings: [[1]] = outer ring, [[2+]] = holes
      if (length(poly) > 1) {
        holes <- poly[-1]          # drop the outer ring
        for (hole in holes) {
          # A hole ring is counter-clockwise; convert to a proper polygon
          island_rings <- c(island_rings, list(sf::st_polygon(list(hole))))
        }
      }
    }
  }
  
  if (length(island_rings) == 0) {
    message("No islands found.")
    return(NULL)
  }
  
  islands_sfc <- sf::st_sfc(island_rings, crs = sf::st_crs(x))
  sf::st_sf(geometry = islands_sfc)
}


#' Extract the outer shoreline of a lake multipolygon
#'
#' Returns the exterior boundary of the lake, ignoring any island holes.
#' The result is a POLYGON (or MULTIPOLYGON) with no holes.
#'
#' @param x An sf object with POLYGON or MULTIPOLYGON geometry
#' @return An sf object with the outer shoreline polygon(s)
#' @export
#' @importFrom sf st_as_sf
#' @importFrom terra vect fillHoles
detect_shoreline <- function(x) {
  stopifnot(inherits(x, "sf") || inherits(x, "sfc"))
  
  x |>
    sf::st_geometry() |> 
    terra::vect() |>
    terra::fillHoles() |>
    # terra::as.lines() |> 
    sf::st_as_sf()
}

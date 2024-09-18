#' Generate depth points for interpolation
#'
#' @param shoreline sf object of lake shoreline.
#' @param islands sf object of lake islands if present.
#' @param point_data sf object of depth points. Must contain a "depth" column.
#' If NULL, then contours must be provided.
#' @param contours sf object of depth contours. If NULL, then points must be
#' provided.
#' @param res numeric resolution of output raster in metres.
#' @param subsample logical subsample points to 10000 if TRUE. Default is TRUE.
#' @inheritParams sf::st_transform
#'
#' @importFrom sf st_cast st_length st_area st_transform st_intersection st_crs
#' st_sample st_as_sf st_buffer st_line_sample st_make_valid st_is_empty
#' @importFrom units set_units drop_units
#' @importFrom dplyr group_by filter mutate select bind_rows sample_n
#'
#' @return sf object of depth points.
#' @export
#'
#' @examples
#' shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
#' package = "bathytools"))
#' point_data <- readRDS(system.file("extdata/depth_points.rds",
#' package = "bathytools"))
#'
#' # Generate depth points with points ----
#' depth_points <- generate_depth_points(shoreline = shoreline,
#' point_data = point_data)
#'
#' # Generate depth points with contours ----
#' contours <- readRDS(system.file("extdata/depth_contours.rds",
#' package = "bathytools"))
#' depth_points <- generate_depth_points(shoreline = shoreline,
#' contours = contours)
#'

generate_depth_points <- function(shoreline, islands = NULL, point_data = NULL,
                                  contours = NULL, res = 2, subsample = TRUE,
                                  crs = NULL) {

  message("Generating depth points... [", format(Sys.time()), "]")

  if (is.null(crs)) {
    crs <- sf::st_crs(shoreline)
  }

  if (!is.null(islands)) {
    # Add handler for islands

  }

  # If shoreline CRS is not equal to crs then transform
  if (sf::st_crs(shoreline) != crs) {
    shoreline <- sf::st_transform(shoreline, crs)
  }

  # If shoreline is a polygon, convert to line
  shoreline_poly <- shoreline
  shoreline_line <- shoreline
  if (sf::st_geometry_type(shoreline) %in% c("POLYGON", "MULTIPOLYGON")) {
    geom <- shoreline |>
      sf::st_geometry() |>
      sf::st_cast("MULTILINESTRING") |>
      sf::st_cast("LINESTRING")
    sf::st_geometry(shoreline_line) <- geom
    sf::st_geometry(shoreline_line) <- "geometry"
  } else if (sf::st_geometry_type(shoreline) %in% c("MULTILINESTRING", "LINESTRING")) {
    shoreline_poly <- shoreline |>
      sf::st_cast("POLYGON")
  }



  # Calculate lake perimeter for shoreline ----
  perim <- shoreline_line |>
    sf::st_length() |>
    units::set_units("m") |>
    units::drop_units()
  nperim <- round(perim / res)
  if (nperim > 1000) {
    nperim <- 500
  }

  # Calculate lake area
  lake_area <- shoreline_poly |>
    sf::st_geometry() |>
    sf::st_area() |>
    units::set_units("m^2") |>
    units::drop_units()



  if (!is.null(contours)) {

    if (!"depth" %in% names(contours)) {
      stop("Contours must have a column named 'depth'")
    }
    # Check if all contour depths are less than or equal to 0
    if (any(contours$depth > 0)) {
      stop("Contour depths must be less than or equal to 0")
    }

    # Transform contours to output CRS
    contours <- contours |>
      sf::st_transform(crs)


    # Compare istance bet
    grp <- contours |>
      # sf::st_make_valid() |>
      dplyr::group_by(depth)
    chk <- sapply(1:nrow(grp), function(i) {
      tryCatch({
        grp[i,] |>
          sf::st_cast("MULTIPOINT") |>
          sf::st_cast("POLYGON")
        TRUE
      }, error = function(e) {
        return(FALSE)
      })
    })
    conts <- grp[chk, ] |>
      sf::st_cast("MULTIPOINT") |>
      sf::st_cast("POLYGON") |>
      dplyr::reframe(area = sf::st_area(geometry)) |>
      # dplyr::summarise(area = sum(area)) |>
      dplyr::mutate(area = units::drop_units(area))
    if (any(duplicated(conts$depth))) {
      conts <- conts |>
        dplyr::group_by(depth) |>
        dplyr::summarise(area = sum(area))
    }

    # Sample points regularly for each depth polygon
    cont_pnts <- lapply(unique(contours$depth), function(d) {

      perim <- contours |>
        dplyr::filter(depth == d) |>
        sf::st_cast("MULTILINESTRING") |>
        sf::st_length() |>
        units::set_units("m") |>
        units::drop_units()
      np <- round(perim / res)

      p <- contours |>
        dplyr::filter(depth == d) |>
        sf::st_cast("MULTILINESTRING") |>
        sf::st_cast("LINESTRING") |>
        sf::st_line_sample(n = 10, type = "regular") |>
        # sf::st_sample(np, type = "regular") |>
        sf::st_as_sf() |>
        sf::st_cast("POINT") |>
        dplyr::mutate(depth = d) |>
        dplyr::rename(geometry = x)

      line_string <- contours |>
        dplyr::filter(depth == d) |>
        sf::st_cast("MULTILINESTRING") |>
        sf::st_cast("LINESTRING")

      p <- lapply(1:nrow(line_string), \(i) {
        lgth <- units::drop_units(sf::st_length(line_string[i, ]))
        np <- ceiling(lgth / res)

        line_string[i, ] |>
          sf::st_line_sample(n = np, type = "regular") |>
          # sf::st_sample(np, type = "regular") |>
          sf::st_as_sf() |>
          sf::st_cast("POINT") |>
          dplyr::mutate(depth = d) |>
          dplyr::rename(geometry = x)
      }) |>
        dplyr::bind_rows()
      p
    }) |>
      dplyr::bind_rows() |>
      sf::st_as_sf() |>
      dplyr::select(depth, geometry)
  }

  # If point_data dataframe is provided, convert to sf object
  if (!is.null(point_data)) {

    if (!"depth" %in% names(point_data)) {
      stop("point_data must have a column named 'depth'")
    }
    if (any(point_data$depth > 0)) {
      stop("Point depths must be less than or equal to 0")
    }

    # If point_data are a dataframe, convert to sf object
    if (is.data.frame(point_data)) {

      if (!all(c("lon", "lat", "depth") %in% names(point_data))) {
        stop("point_data must have columns 'lon', 'lat', and 'depth'")
      }

      pts <- sf::st_as_sf(point_data, coords = c("lon", "lat"),
                          crs = 4326) |>
        dplyr::mutate(depth = point_data$depth) |>
        dplyr::select(depth, geometry)
    } else if (is(sf::st_geometry(point_data), "POINT")) {
      pts <- point_data |>
        dplyr::select(depth, geometry)
    }
    pts <- sf::st_transform(pts, crs)
  }

  if (!is.null(contours) & !is.null(point_data)) {
    point_data <- dplyr::bind_rows(cont_pnts, pts)
  } else if (!is.null(contours)) {
    point_data <- cont_pnts
  } else if (!is.null(point_data)) {
    point_data <- pts
  }

  # Generate an outline of the lake as points for interpolation ----
  pnt_outline <- shoreline_line |>
    sf::st_sample(nperim, type = "regular") |>
    sf::st_as_sf() |>
    sf::st_cast("POINT") |>
    dplyr::rename(geometry = x) |>
    dplyr::mutate(depth = 0) |>
    dplyr::select(depth, geometry)

  # Remove empty geometry
  if (any(sf::st_is_empty(pnt_outline))) {
    pnt_outline <- pnt_outline |>
      dplyr::filter(!sf::st_is_empty(geometry))
  }

  # Generate points around lake perimeter for interpolation ----
  pnt_inset <- shoreline_poly |>
    sf::st_buffer(0.1) |>
    sf::st_cast("MULTILINESTRING") |>
    sf::st_sample(nperim, type = "regular") |>
    sf::st_as_sf() |>
    sf::st_cast("POINT") |>
    dplyr::rename(geometry = x) |>
    dplyr::mutate(depth = -0.1) |>
    dplyr::select(depth, geometry)

  if (any(sf::st_is_empty(pnt_inset))) {
    pnt_inset <- pnt_inset |>
      dplyr::filter(!sf::st_is_empty(geometry))
  }

  if (subsample) {
    if (nrow(point_data) > 10000) {
      message("Warning: large number of points for interpolation (", nrow(point_data), ")")
      # Subsample to 10000
      point_data <- point_data |>
        dplyr::sample_n(10000, replace = FALSE)
    }
  }

  # Subset point_data to those within the lake
  # sub1 <- point_data[sf::st_within(point_data, shoreline_poly, sparse = FALSE), ]
  suppressWarnings({
    sub <- point_data |>
      sf::st_intersection(shoreline_poly) |>
      dplyr::select(depth)
  })

  # Combine in-lake point_data & outline
  all <- dplyr::bind_rows(pnt_outline, pnt_inset, sub) |>
    dplyr::mutate(depth = depth)


  if (nrow(all) > 15000) {
    message("Warning: large number of points for interpolation (", nrow(all), ")")
    if (subsample) {
      message("Subsampling to 10000 points")
      # Subsample to 10000
      all <- all |>
        dplyr::sample_n(10000)
    }
  }
  message("Finished! [", format(Sys.time()), "]")
  return(all)
}

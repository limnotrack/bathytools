#' Extend a hypsograph from a DEM
#' @description This function calculates a hypsograph from a bathymetry merged
#'  with a DEM, optionally extending it.
#' to include an area above the lake surface. It can also calculate the area of
#' the lake at different depths, and optionally mask the DEM to a shoreline.
#' @inheritParams merge_bathy_dem
#' @inheritParams get_depths
#' @param shoreline sf object. The shoreline polygon to mask the DEM to.
#' @param dem_bath SpatRaster object. The DEM with bathymetry data.
#' @param lake_elev numeric. The elevation of the lake surface. If NULL,
#' the function will calculate the lake surface elevation from the DEM.
#' @param ext_elev numeric. The elevation above the lake surface to extend the
#' hypsograph. Default is 0, which means no extension.
#' @return data.frame with columns 'elev', 'area', and 'depth'. The 'elev'
#' column contains the elevation, the 'area' column contains the area of the
#' lake at that elevation, and the 'depth' column contains the depth of the lake
#' at that elevation.
#' @importFrom terra minmax res ifel mask expanse
#' @importFrom dplyr arrange mutate bind_rows
#' @importFrom sf st_as_sf st_cast
#' @importFrom tmap tm_shape tm_raster tm_borders
#' @importFrom utils data
#' @export

dem_to_hypsograph <- function(shoreline = NULL, dem_bath, lake_elev = NULL,
                              lake_depth = NULL, 
                              depths = NULL, ext_elev = 0) {
  
  # Mask to shoreline if provided
  if (!is.null(shoreline)) {
    in_lake <- terra::mask(dem_bath, shoreline)
    out_lake <- terra::mask(dem_bath, shoreline, inverse = TRUE)
  } else {
    in_lake <- dem_bath
    out_lake <- NULL
  }
  
  # Basic raster properties
  cell_area <- prod(terra::res(dem_bath))
  mm <- terra::minmax(dem_bath)
  
  min_lake_elev <- terra::minmax(in_lake)[1]
  max_lake_elev <- terra::minmax(in_lake)[2]
  
  if (is.null(lake_elev)) {
    lake_elev <- max_lake_elev
  }
  
  depth_range <- lake_elev - min_lake_elev
  if (is.null(lake_depth)) {
    lake_depth <- max_lake_elev - min_lake_elev
    message("Lake depth not provided, using maximum lake depth: ", lake_depth)
  } else {
    lake_depth <- abs(lake_depth)
  }
  
  # Depth interval defaults
  if (is.null(depths)) {
    
    if (!lake_depth %in% model_layer_structure$zi) {
      mod_layers <- model_layer_structure |> 
        dplyr::bind_rows(
          data.frame(zi = lake_depth, z = lake_depth)
        ) |> 
        dplyr::arrange(zi)
    }
    depths <- mod_layers |> 
      dplyr::filter(z <= lake_depth) |> 
      dplyr::mutate(depths = lake_elev - zi) |> 
      dplyr::pull(depths)
  }
  
  # Core hypsograph below lake surface
  hyps <- data.frame(
    elev = depths,
    area = sapply(depths, function(e) sum(in_lake[] <= e, na.rm = TRUE) * cell_area)
  )
  hyps$depth <- -round(lake_elev - hyps$elev, 3)
  hyps <- dplyr::arrange(hyps, dplyr::desc(elev))
  # plot(hyps$area, hyps$elev, main = "Hypsograph Below Lake Surface")
  
  # Optional extension above lake surface
  if (ext_elev > 0) {
    top_elev <- lake_elev + ext_elev
    if (top_elev > mm[2]) {
      message("Top elevation exceeds DEM max — truncated.")
      top_elev <- mm[2]
    }
    
    poly_max <- extract_ext_elev_polygon(dem_bath = dem_bath, 
                                         lake_elev = lake_elev, 
                                         ext_elev = ext_elev, 
                                         shoreline = shoreline)
    out_lake_ext <- terra::mask(dem_bath, poly_max)
    terra::plot(out_lake_ext, main = "Extended Lake Area")
    
    ext_depths <- model_layer_structure |> 
      dplyr::mutate(depths = lake_elev + zi) |> 
      dplyr::filter(depths < top_elev, zi > 0) |> 
      dplyr::pull(depths)
    
    ext_hyps <- data.frame(
      elev = ext_depths,
      area = sapply(ext_depths, function(e) sum(out_lake_ext[] <= e,
                                                na.rm = TRUE) * cell_area),
      depth = ext_depths - lake_elev
    )
    # plot(hyps$area, hyps$elev, main = "Hypsograph Below Lake Surface")
    # points(ext_hyps$area, ext_hyps$elev, col = "red", pch = 19)
    # # ext_hyps$depth <- round(lake_elev - ext_hyps$elev, 3)
    # 
    hyps <- dplyr::bind_rows(hyps, ext_hyps) |>
      dplyr::arrange(dplyr::desc(elev)) |>
      dplyr::select(elev, depth, area)
    # plot(hyps$area, hyps$elev, main = "Hypsograph with Extension")
    # abline(h = lake_elev, col = "red", lwd = 2)
  }
  
  return(hyps)
}

#' Extract the largest polygon containing the lake from an extended elevation area
#' @description This function extracts the largest polygon that intersects the
#' shoreline polygon from a raster representing an extended elevation area above
#' the lake surface.
#' @param dem_bath SpatRaster object. The DEM with bathymetry data.
#' @param lake_elev numeric. The elevation of the lake surface.
#' @param ext_elev numeric. The elevation above the lake surface to extend the
#' hypsograph.
#' @param shoreline sf object. The shoreline polygon to intersect with.
#' @return sf object. The largest polygon containing the lake from the extended
#' elevation area.
#' @importFrom terra ifel as.polygons
#' @importFrom sf st_as_sf st_cast st_intersects st_transform
#' @importFrom dplyr filter mutate bind_rows
#' @noRd

extract_ext_elev_polygon <- function(dem_bath, lake_elev, ext_elev, shoreline) {
  # Create raster mask of cells <= lake_elev + ext_elev
  mask_rast <- terra::ifel(dem_bath <= (lake_elev + ext_elev), 1, NA)
  
  # Convert to polygons
  polys <- terra::as.polygons(mask_rast, dissolve = TRUE) |>
    sf::st_as_sf() |>
    sf::st_cast("POLYGON")
  
  # Ensure shoreline is in same CRS
  if (sf::st_crs(polys) != sf::st_crs(shoreline)) {
    shoreline <- sf::st_transform(shoreline, sf::st_crs(polys))
  }
  
  # Filter polygons that intersect the shoreline polygon
  intersects <- sf::st_intersects(polys, shoreline, sparse = FALSE)[, 1]
  
  if (!any(intersects)) {
    stop("No polygon containing lake found in extended elevation area.")
  }
  
  # Among polygons intersecting shoreline, pick the largest by area
  polys_sub <- polys[intersects, ]
  areas <- sf::st_area(polys_sub)
  largest_poly <- polys_sub[which.max(areas), ]
  
  return(largest_poly)
}


# 
# 
# dem_to_hypsograph <- function(shoreline = NULL, dem_bath, lake_elev = NULL,
#                               depths = NULL, ext_elev = 0) {
#   
#   utils::data("model_layer_structure", package = "bathytools",
#               envir = environment())
#   
#   if (!is.null(shoreline)) {
#     in_lake <- terra::mask(dem_bath, shoreline)
#     out_lake <- terra::mask(dem_bath, shoreline, inverse = TRUE)
#   } else {
#     in_lake <- dem_bath
#   }
#   
#   rast_res <- raster::res(dem_bath)
#   cell_area <- rast_res[1] * rast_res[2]
#   mm <- terra::minmax(dem_bath)
#   max_elev <- round(mm[2], 2)
#   min_elev <- round(mm[1], 2)
#   min_lake_elev <- round(terra::minmax(in_lake)[1], 2)
#   max_lake_elev <- round(terra::minmax(in_lake)[2], 2)
#   if (is.null(lake_elev)) {
#     lake_elev <- max_lake_elev
#   }
#   # if (lake_elev != max_lake_elev) {
#   #   message("shoreline elevation (", lake_elev, "m) is not equal to shoreline elevation from raster (", max_lake_elev, "m)")
#   # }
#   depth <- lake_elev - min_lake_elev
#   
#   # Set dem_bath to NA below lake_elev
#   in_lake <- terra::ifel(in_lake < min_lake_elev, NA, in_lake)
#   in_lake <- terra::ifel(in_lake > lake_elev, lake_elev, in_lake)
#   # terra::plot(in_lake)
#   # Add intervals if range is less than 10m and greater than 100m ----
#   if (depth <= 5) {
#     intv <- 0.2
#   } else if (depth <= 50) {
#     intv <- 0.5
#   } else if (depth <= 100) {
#     intv <- 1
#   } else {
#     intv <- 5
#   }
#   
#   top_elev <- lake_elev + ext_elev
#   if (top_elev > max_elev) {
#     message("Top elevation is greater than maximum elevation of the DEM.\nSetting top elevation to maximum elevation.")
#     top_elev <- max_elev
#     ext_elev <- 0
#   }
#   
#   if (is.null(depths)) {
#     depths <- model_layer_structure |> 
#       dplyr::filter(z <= depth) |> 
#       dplyr::mutate(depths = lake_elev - zi) |> 
#       dplyr::pull(depths)
#     # depths <- seq(min_lake_elev, lake_elev, by = intv)
#     # depths <- get_depths(bathy_raster = dem_bath, surface = lake_elev, depths = intv)
#   }
#   # if (!(lake_elev %in% depths)) {
#   #   depths <- c(depths, lake_elev)
#   # }
#   depths <- depths[order(depths)]
#   
#   # Calculate hypsograph ----
#   l_area <- sum(!is.na(dem_bath[])) * (cell_area)
#   # chk <- dem_bath <= -210
#   # l_area <- terra::expanse(dem_bath)
#   # l_area <- l_area[2]
#   hyps <- data.frame(elev = depths,
#                      area = sapply(depths, \(d) sum(in_lake[] <= d, na.rm = TRUE)),
#                      depth = lake_elev - depths) |> 
#     dplyr::mutate(area = area * (cell_area),
#                   depth = -1 * depth
#     ) |> 
#     dplyr::arrange(dplyr::desc(elev))
#   hyps
#   plot(hyps$area, hyps$elev)
#   
#   if (ext_elev > 0) {
#     # Extract the extended elevation polygon
#     poly_max <- extract_ext_elev_polygon(dem_bath = dem_bath, 
#                                          shoreline = shoreline,
#                                          lake_elev, ext_elev)
#     tm_shape(dem_bath) +
#       tm_raster() +
#       tm_shape(poly_max) +
#       tm_borders(col = "red", lwd = 2) +
#       tm_shape(shoreline) +
#       tm_borders(col = "blue", lwd = 2)
#     
#     out_lake_ext <- terra::mask(dem_bath, poly_max)
#     
#     terra::plot(out_lake_ext)
#     
#     
#     depths <- seq(lake_elev, top_elev, by = intv)
#     depths <- depths[-1]
#     if (!(top_elev %in% depths)) {
#       depths <- c(depths, top_elev)
#     }
#     depths <- depths[order(depths)]
#     hyps2 <- data.frame(elev = depths,
#                         area = sapply(depths, \(d) sum(out_lake_ext[] <= d, na.rm = TRUE)),
#                         depth = lake_elev - depths) |> 
#       dplyr::mutate(area = area * (cell_area),
#                     depth = -1 * depth,
#                     # norm_depth = NA, 
#                     # norm_area = NA
#       ) |> 
#       dplyr::arrange(dplyr::desc(elev))
#     plot(hyps2$area, hyps2$elev)
#     hyps <- dplyr::bind_rows(hyps, hyps2) |> 
#       dplyr::arrange(dplyr::desc(elev))
#   }
#   hyps <- hyps |> 
#     dplyr::mutate(depth = round(depth, 3))
#   
#   plot(hyps$area, hyps$elev)
#   abline(h = lake_elev, col = "red")
#   
#   return(hyps)
# }
# 
# extract_ext_elev_polygon <- function(dem_bath, lake_elev, ext_elev) {
#   # Extract largest polygon connected to the shoreline for calculating extended hypsograph
#   test <- terra::clamp(dem_bath, upper = (lake_elev + ext_elev), values = FALSE)
#   test <- terra::ifel(dem_bath > (lake_elev + ext_elev), NA, lake_elev + ext_elev)
#   # terra::plot(test)
#   poly <- terra::as.polygons(test) |> 
#     sf::st_as_sf() |> 
#     sf::st_cast("POLYGON") 
#   
#   p_area <- poly |> 
#     sf::st_area()
#   poly_max <- poly[p_area == max(p_area), ]
#   return(poly_max)
#   
# }

#' Extract the largest polygon connected to the lake for calculating extended hypsograph
#' @param dem_bath SpatRaster object with the bathymetry data.
#' @param shoreline sf object with the shoreline of the lake.
#' @param lake_elev numeric. The elevation of the lake.
#' @param ext_elev numeric. The additional elevation to extend the lake surface.
#' @return sf object with the largest polygon connected to the lake.
#' @noRd
#' @importFrom terra clamp ifel as.polygons disagg subset relate vect project 
#' crs
#' @importFrom sf st_as_sf st_cast st_intersects
#' @importFrom dplyr filter
extract_ext_elev_polygon <- function(dem_bath, shoreline, lake_elev, ext_elev) {
  # Extract largest polygon connected to the lake for calculating extended 
  # hypsograph
  # --- build a mask (1 = keep, NA = drop) ---
  mask <- dem_bath <= (lake_elev + ext_elev)      # logical SpatRaster
  mask <- terra::ifel(mask, 1, NA)                # 1 / NA
  
  # --- convert to polygons, merging contiguous cells ---
  ext_poly <- terra::as.polygons(mask, dissolve = TRUE) |> 
    terra::disagg()
  
  shoreline_vect <- shoreline |> 
    terra::vect() |> 
    terra::project(terra::crs(dem_bath))
  
  # plot(shoreline, add = TRUE)
  # ext_poly <- terra::as.polygons(dem_cond)
  # terra::plot(ext_poly)
  # terra::plot(shoreline_vect, add = TRUE)
  # 
  # Check intersection
  rel <- terra::relate(ext_poly, shoreline_vect, relation = "intersects")

  # Subset polygons
  poly <- ext_poly |> 
    terra::subset(rel[, 1])
  
  # terra::plot(poly, add = TRUE, border = "red", lwd = 2)
  
  return(poly)
}


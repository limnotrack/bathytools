extract_ext_elev_polygon <- function(dem_bath, shoreline, lake_elev, ext_elev) {
  # Extract largest polygon connected to the lake for calculating extended hypsograph
  dem_cond <- terra::clamp(dem_bath, upper = (lake_elev + ext_elev), values = FALSE)
  dem_cond <- terra::ifel(dem_bath > (lake_elev + ext_elev), NA, lake_elev + ext_elev)
  # terra::plot(dem_cond)
  plot(shoreline, add = TRUE)
  ext_poly <- terra::as.polygons(dem_cond) |> 
    sf::st_as_sf() |> 
    sf::st_cast("POLYGON") 
  
  # Find which polygon contains the shoreline
  poly <- ext_poly[sf::st_intersects(ext_poly, shoreline, sparse = FALSE), ]
  # terra::plot(poly)
  
  return(poly)
}

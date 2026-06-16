#' Draw a elevation map with bathymetry and surrounding topography
#'
#' @param dem_bath SpatRaster object with merged DEM and bathymetry data. Can be
#' generated using the \link{merge_bathy_dem} function.
#' @param lake_elev numeric, elevation of the lake surface. Can be extracted using
#' the \link{get_lake_surface_elevation} function.
#' @param ext_elev numeric, additional elevation to extend the lake surface.
#' This allows you to generate a flood map. Default is 0.
#'
#' @inherit tmap::tm_shape return
#'
#' @importFrom terra minmax
#' @importFrom tmap tm_shape tm_raster
#' @export
#'

tm_dem_bath <- function(dem_bath, lake_elev, ext_elev = 0) {

  if (ext_elev == 0) {
    ext_elev <- 0.1 # prevents duplication
  }

  mm <- terra::minmax(dem_bath)
  lake_cuts <- round(seq(from = mm[1], to = lake_elev, length.out = 5), 1)
  elev_cuts <- round(seq(from = (lake_elev + ext_elev), to = mm[2], length.out = 5), 1)
  elev_cuts <- elev_cuts[-1]
  cuts <- c(lake_cuts, elev_cuts) #set breaks
  flood_cuts <- c(lake_elev + ext_elev)
  cuts <- c(lake_cuts, flood_cuts, elev_cuts) #set breaks

  lake_surf <- (lake_elev - mm[1]) / (mm[2] - mm[1])
  flood1 <- (lake_elev - mm[1] + 0.1) / (mm[2] - mm[1])
  flood2 <- (lake_elev - mm[1] + ext_elev) / (mm[2] - mm[1])
  terr_values <- seq(flood2+0.01, 1, length.out = 6)
  values <- c(0, lake_surf, flood1, terr_values)
  breaks <- round(values * (mm[2] - mm[1]) + mm[1], 1)
  breaks <- c(floor(mm[1]), lake_elev, lake_elev + ext_elev,
              round(terr_values[-1] * (mm[2] - mm[1]) + mm[1], 1))
  breaks[length(breaks)] <- mm[2]
  colours <- c(c("#084594", "white", "#C7E9C0", "#41AB5D", "#E6E600",
                 "#EBB25E", "#F0C9C0"))

  tmap::tm_shape(dem_bath) +
    tmap::tm_raster(col.scale = tmap::tm_scale_continuous(values = colours, 
                                                          ticks = breaks))
}

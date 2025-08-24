#' Estimate bathymetry from a shoreline polygon using distance transformation.
#'
#' This function estimates the bathymetry of a lake from a shoreline polygon using
#' a distance transformation. The function calculates the distance from the
#' shoreline polygon to the edge of the lake and uses this distance to estimate
#' the bathymetry. The maximum depth of the lake is provided as an input parameter.
#' This is the formula used in the Hollister and Milstead (2010) paper.
#'
#' @inheritParams bathy_to_hypso
#' @param max_depth numeric. The maximum depth of the lake.
#' @param hypsograph data.frame. A hypsograph data frame with columns 'depth' and
#' 'area'. If provided, the function will use the hypsograph to estimate the
#' bathymetry. If not provided, the function will use a linear scaling based on
#' the maximum depth.
#' @inheritParams generate_depth_points
#'
#' @return SpatRaster object with the estimated bathymetry.
#' @export
#'
#' @references Hollister, J. W., W.B. Milstead (2010). Using GIS to Estimate
#'             Lake Volume from Limited Data. Lake and Reservoir Management.
#'             26(3)194-199.
#'             \doi{10.1080/07438141.2010.504321}
#' @examples
#' shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
#' package = "bathytools"))
#' max_depth <- 81.6
#' bathy <- estimate_bathymetry(shoreline = shoreline, max_depth = max_depth)
#' terra::plot(bathy)

estimate_bathymetry <- function(shoreline, hypsograph = NULL, max_depth, 
                                res = 10) {
  
  lake_vect <- terra::vect(shoreline)
  lake_extent <- terra::ext(lake_vect)
  # Create a blank raster based on the extent and resolution
  distance_raster <- terra::rast(ext = lake_extent, resolution = res,
                                 crs = terra::crs(lake_vect))
  terra::values(distance_raster) <- 1
  
  # Rasterize the lake polygon
  # Assign a value of 1 to the lake polygon area
  lake_mask <- terra::mask(distance_raster, lake_vect, updatevalue = NA,
                           inverse = TRUE)
  
  # Calculate the distance from the lake polygon
  # Distance will be measured in the same units as the coordinate system of the lake_polygon
  distance_raster <- terra::distance(lake_mask) |> 
    terra::mask(lake_vect)
  
  if (!is.null(hypsograph)) {
    # 1) Clean & normalize hypsograph
    hypsograph <- hypsograph |> 
      dplyr::filter(depth <= 0) |>  # Keep only depths below surface
      dplyr::arrange(-depth)
    A0 <- max(hypsograph[["area"]], na.rm = TRUE)
    hypsograph <- hypsograph |> 
      dplyr::mutate(area_frac = area / A0) # Normalize area
    
    
    # 2) Collect lake-cell distances
    dvals <- terra::values(distance_raster, mat = FALSE, na.rm = TRUE)
    if (!length(dvals)) stop("distance_raster has no in-lake (non-NA) cells.")
    
    # 3) For each depth, find distance cutoff whose enclosed area matches hyps area
    # If F(t) = P(distance <= t), then area beyond t is 1 - F(t).
    # We want 1 - F(t_d) = area_frac(depth d)  =>  F(t_d) = 1 - area_frac
    probs <- 1 - hypsograph$area_frac
    
    # Quantiles of the empirical distance distribution
    # type=7 is default; robust choices are 1 or 2 if you want stepwise behavior
    cuts <- stats::quantile(dvals, probs = probs, na.rm = TRUE, type = 7,
                            names = FALSE)
    
    # 4) Build a monotone map distance -> depth using these (distance_cutoff, depth) pairs
    # (distance increases from shore to interior; depth should increase in magnitude)
    map_fun <- stats::approxfun(x = cuts, y = hypsograph[["depth"]], rule = 2, 
                                ties = "ordered")
    
    # 5) Apply to the raster
    bathymetry_raster <- terra::app(distance_raster,
                                    fun = function(x) map_fun(x))
    terra::plot(bathymetry_raster)
    
    
  } else {
    max_depth <- ifelse(max_depth > 0, -max_depth, max_depth)
    adj_depth <- round(0.01 * max_depth, 2)
    max_distance <- max(terra::values(distance_raster), na.rm = TRUE)
    # Apply a function to distance raster
    bathymetry_raster <- terra::app(distance_raster, function(x) {
      depth = (x * max_depth) / max_distance
    })
  }
  
  bathymetry_raster <- terra::mask(bathymetry_raster, lake_vect)
  # bathymetry_raster <- -bathymetry_raster
  names(bathymetry_raster) <- "depth"
  return(bathymetry_raster)
}

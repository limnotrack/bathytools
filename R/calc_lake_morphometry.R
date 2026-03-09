#' Calculate comprehensive morphometric statistics from a lake bathymetry raster
#'
#' @description
#' Computes a comprehensive set of lake morphometric statistics from a
#' bathymetric raster. Metrics span six categories: basic geometry, depth
#' statistics, morphometric ratios, shoreline shape, basin slope, and
#' hypsographic curve parameters. Together these characterise the 3-D form of
#' the lake basin and are suitable as features for predictive modelling or
#' comparative limnological analysis.
#'
#' **1. Basic Geometry**
#'
#' - `surface_area_m2` / `surface_area_km2`: Planimetric area of the water
#'   surface, computed as the count of non-NA raster cells multiplied by the
#'   mean cell area. Units: m² and km².
#'
#' - `volume_m3`: Total lake volume computed by summing the product of depth
#'   and cell area across all cells:
#'   \deqn{V = \sum_{i=1}^{n} z_i \cdot a}
#'   where \eqn{z_i} is the depth of cell \eqn{i} and \eqn{a} is the cell area
#'   (m²).
#'
#' - `perimeter_m`: Length of the lake shoreline polygon in metres, derived
#'   either from the supplied `shoreline_vect` or vectorised from the raster
#'   extent.
#'
#' **2. Depth Statistics**
#'
#' All depth statistics are computed from the vector of per-cell depth values.
#'
#' - `z_max`: Maximum depth (m). The deepest recorded point in the basin.
#'
#' - `z_mean`: Arithmetic mean depth (m):
#'   \deqn{\bar{z} = \frac{1}{n} \sum_{i=1}^{n} z_i}
#'   Also computed independently as \eqn{V / A} (`z_mean_from_volume`) as a
#'   consistency check; both values should be nearly identical.
#'
#' - `z_median`: Median depth (m). More robust than the mean when the depth
#'   distribution is strongly skewed.
#'
#' - `z_sd`: Standard deviation of depth (m):
#'   \deqn{s_z = \sqrt{\frac{1}{n-1} \sum_{i=1}^{n} (z_i - \bar{z})^2}}
#'
#' - `z_cv`: Coefficient of variation of depth (unitless). Normalises
#'   variability by the mean, enabling comparison across lakes of different
#'   sizes:
#'   \deqn{CV_z = \frac{s_z}{\bar{z}}}
#'
#' - `z_skewness`: Pearson's moment coefficient of skewness (type 2,
#'   bias-corrected; via `e1071::skewness`). Positive values indicate that most
#'   of the lake area is shallow with a long tail of deeper water (common in
#'   glacially scoured basins); negative values indicate a flat bottom with
#'   steep walls.
#'   \deqn{g_1 = \frac{n}{(n-1)(n-2)} \sum_{i=1}^{n}
#'         \left(\frac{z_i - \bar{z}}{s_z}\right)^3}
#'
#' - `z_kurtosis`: Excess kurtosis of the depth distribution (type 2,
#'   bias-corrected; via `e1071::kurtosis`). High positive kurtosis indicates
#'   a peaked distribution (many cells near a single dominant depth); values
#'   near zero indicate a flat, uniform distribution.
#'
#' - `z_p10`, `z_p25`, `z_p50`, `z_p75`, `z_p90`: Depth percentiles (m) at
#'   the 10th, 25th, 50th, 75th, and 90th quantiles of the depth distribution.
#'   Together these provide a non-parametric summary of the hypsographic curve
#'   without requiring a functional fit.
#'
#' **3. Morphometric Ratios**
#'
#' - `volume_development` (\eqn{V_d}): Compares the actual basin volume to the
#'   volume of a cone with the same surface area and maximum depth. A cone has
#'   \eqn{V_d = 1}; values < 1 indicate a concave (bowl-shaped) basin where
#'   volume is concentrated near the surface; values > 1 indicate a convex
#'   (flat-bottomed) basin. Introduced by Hutchinson (1957):
#'   \deqn{V_d = \frac{3\,\bar{z}}{z_{max}}}
#'
#' - `relative_depth` (\eqn{Z_r}): Expresses maximum depth as a percentage of
#'   the mean diameter of the lake surface, allowing depth to be compared across
#'   lakes of very different sizes (Hutchinson 1957):
#'   \deqn{Z_r = \frac{50\, z_{max} \sqrt{\pi}}{\sqrt{A}}}
#'   where \eqn{A} is surface area in m². High values indicate a deep,
#'   small-surface lake; low values a shallow, large-surface lake. Strongly
#'   linked to stratification stability.
#'
#' - `dynamic_ratio` (\eqn{DR}): A predictor of thermal stratification and wind
#'   mixing (Håkanson 1982). High values indicate lakes prone to full mixing;
#'   low values indicate stable stratification:
#'   \deqn{DR = \frac{\sqrt{A_{km^2}}}{\bar{z}}}
#'   where \eqn{A_{km^2}} is surface area in km² and \eqn{\bar{z}} is mean
#'   depth in metres. Lakes with \eqn{DR > \approx 0.8} are generally
#'   considered unstratified.
#'
#' **4. Shoreline & Shape Metrics**
#'
#' - `shoreline_development` (\eqn{D_L}): Ratio of the actual shoreline length
#'   to the circumference of a circle enclosing the same area. A perfect circle
#'   gives \eqn{D_L = 1}; larger values indicate greater shoreline irregularity
#'   and thus greater littoral habitat complexity (Hutchinson 1957):
#'   \deqn{D_L = \frac{L}{2\sqrt{\pi A}}}
#'   where \eqn{L} is perimeter (m) and \eqn{A} is surface area (m²).
#'
#' - `circularity`: The isoperimetric quotient, an alternative shape index
#'   bounded on \\[0, 1\\] where 1 is a perfect circle. Mathematically the
#'   reciprocal of \eqn{D_L^2}:
#'   \deqn{C = \frac{4\pi A}{L^2}}
#'
#' - `elongation_ratio`: Ratio of the longest to the shortest axis of the
#'   axis-aligned bounding box of the shoreline polygon. Values near 1 indicate
#'   a compact, equidimensional lake; large values indicate an elongated lake,
#'   which affects fetch, internal seiching, and circulation patterns.
#'   \deqn{E = \frac{\max(\Delta x,\, \Delta y)}{\min(\Delta x,\, \Delta y)}}
#'
#' **5. Basin Slope Metrics**
#'
#' Slope is computed cell-by-cell from the bathymetric raster using
#' `terra::terrain(v = "slope", unit = "degrees")`, which applies a Horn
#' (1981) finite-difference gradient estimator to the eight-cell neighbourhood.
#'
#' - `slope_mean`: Mean basin slope (degrees). Steeper mean slopes are
#'   associated with younger, tectonically active, or glacially carved basins
#'   and with greater sediment focusing toward the profundal zone.
#'
#' - `slope_median`: Median basin slope (degrees). Less sensitive than the mean
#'   to extreme cliff or wall cells at basin margins.
#'
#' - `slope_sd`: Standard deviation of slope (degrees). High values indicate
#'   heterogeneous relief with both very flat and very steep zones.
#'
#' - `slope_skewness`: Skewness of the slope distribution. Positive skew is
#'   typical (most cells are gently sloping, with a tail of steep wall cells);
#'   near-zero skew suggests a uniformly graded basin.
#'
#' - `pct_slope_gentle` / `pct_slope_moderate` / `pct_slope_steep`: Percentage
#'   of total lake area falling into three ecologically relevant slope classes:
#'   gentle (< 5°), moderate (5–15°), and steep (≥ 15°). These classes broadly
#'   correspond to littoral sediment stability zones used in habitat mapping
#'   (Håkanson 1977).
#'
#' **6. Hypsographic Curve Parameters**
#'
#' The hypsographic (or hypsometric) curve describes the cumulative planimetric
#' area of the lake at or above each depth, and is the fundamental 2-D
#' descriptor of basin shape.
#'
#' - `hyps_b_exponent` (\eqn{b}): Exponent of the power-law model fitted to
#'   the normalised hypsographic curve using nonlinear least squares
#'   (`stats::nls`):
#'   \deqn{\hat{A}(z) = A_0 \left(1 - \frac{z}{z_{max}}\right)^b}
#'   where \eqn{\hat{A}(z)} is the area at depth \eqn{z} and \eqn{A_0} is
#'   total surface area. The exponent \eqn{b} summarises basin concavity in a
#'   single parameter: \eqn{b < 1} = concave bowl (volume concentrated near
#'   surface), \eqn{b = 1} = conical basin, \eqn{b > 1} = convex or flat-
#'   bottomed basin. Equivalent to the \eqn{V_d} ratio but estimated
#'   continuously across all depths.
#'
#' - `hyps_r_squared`: Coefficient of determination (R²) for the power-law
#'   fit. Values close to 1 indicate the power model adequately describes the
#'   basin shape; lower values suggest a more complex or irregular hypsograph
#'   that may require a higher-order model.
#'
#' - `hyps_auc_normalised`: Area under the normalised hypsographic curve,
#'   computed by trapezoidal integration over depth values scaled to \\[0, 1\\].
#'   Equivalent to \eqn{\bar{z} / z_{max}} and therefore closely related to
#'   \eqn{V_d / 3}. Values near 1 indicate a flat-bottomed basin; values near
#'   0.5 indicate a conical basin; values near 0 indicate a very concave basin.
#'
#' **7. Depth Zone Areas**
#'
#' - `pct_area_0_2m`, `pct_area_0_5m`, `pct_area_0_10m`: Percentage of total
#'   lake surface area shallower than 2 m, 5 m, and 10 m respectively. These
#'   thresholds correspond approximately to the macrophyte colonisation limit,
#'   the euphotic zone in turbid lakes, and the epilimnion depth in small
#'   temperate lakes.
#'
#' - `pct_area_gt10m`: Percentage of area deeper than 10 m (profundal zone
#'   proxy).
#'
#' - `vol_at_z25pct`, `vol_at_z50pct`, `vol_at_z75pct`: Cumulative water
#'   volume (m³) contained in the shallowest 25%, 50%, and 75% of depth
#'   values respectively. These complement the depth percentiles by weighting
#'   the depth distribution by cell area and reveal how unevenly volume is
#'   distributed through the water column.
#'
#' @references
#' Håkanson, L. (1977). The influence of wind, fetch, and water depth on the
#'   distribution of sediments in Lake Vanern, Sweden. *Canadian Journal of
#'   Earth Sciences*, 14(3), 397–412.
#'
#' Håkanson, L. (1981). *A Manual of Lake Morphometry*. Springer-Verlag, Berlin.
#'
#' Håkanson, L. (1982). Lake bottom dynamics and morphometry: The dynamic
#'   ratio. *Water Resources Research*, 18(5), 1444–1450.
#'
#' Horn, B. K. P. (1981). Hill shading and the reflectance map. *Proceedings
#'   of the IEEE*, 69(1), 14–47.
#'
#' Hutchinson, G. E. (1957). *A Treatise on Limnology. Vol. 1: Geography,
#'   Physics and Chemistry*. Wiley, New York.#'
#' @param bathy_raster A terra::SpatRaster of lake bathymetry. Depth values should
#'   be positive (deeper = larger value) or negative (deeper = more negative).
#'   NA values are treated as outside the lake boundary.
#' @param water_surface_elev Numeric. Elevation of the water surface. Used to
#'   compute depth if raster contains elevation rather than depth. If NULL,
#'   raster values are used directly as depths.
#' @param shoreline Optional terra::SpatVector or sf object of the lake
#'   shoreline polygon. If NULL, the boundary is derived from non-NA raster cells.
#' @param depth_positive Logical. If TRUE (default), deeper water = larger values.
#'   Set FALSE if your raster uses negative values for depth.
#'
#' @return A named list of morphometric statistics
#' 
#' @importFrom cli cli_abort
#'
#' @examples
#' \dontrun{
#'   bathy <- terra::rast("my_lake.tif")
#'   metrics <- calc_lake_morphometry(bathy)
#' }
calc_lake_morphometry <- function(bathy_raster,
                                  water_surface_elev = NULL,
                                  shoreline = NULL,
                                  depth_positive = FALSE) {
  
  # --- 0. Input checks -------------------------------------------------------
  if (!inherits(bathy_raster, "SpatRaster")) {
    cli::cli_abort("bathy_raster must be a terra::SpatRaster object.")
  }
  
  # --- 1. Prepare depth raster -----------------------------------------------
  # Ensure we are working with positive depth values (deeper = larger number)
  depth_rast <- bathy_raster
  
  if (!is.null(water_surface_elev)) {
    # Convert elevation raster to depth
    depth_rast <- water_surface_elev - bathy_raster
  } else if (!depth_positive) {
    # Flip sign so deeper = positive
    depth_rast <- -1 * bathy_raster
  }
  
  # Set negative depths (above water surface / outside lake) to NA
  depth_rast[depth_rast < 0] <- NA
  
  # Extract non-NA depth values as a numeric vector
  depth_vals <- terra::values(depth_rast, na.rm = TRUE) |> as.numeric()
  
  if (length(depth_vals) == 0) {
    stop("No valid (non-NA) depth values found in raster.")
  }
  
  # --- 2. Basic raster properties --------------------------------------------
  cell_area_m2 <- terra::cellSize(depth_rast, unit = "m") |>
    terra::values(na.rm = TRUE) |>
    mean()                          # mean cell area (handles projected CRS)
  
  n_cells      <- length(depth_vals)
  total_area   <- cell_area_m2 * n_cells  # surface area in m²
  
  # --- 3. Depth statistics ---------------------------------------------------
  z_max    <- max(depth_vals)
  z_mean   <- mean(depth_vals)
  z_median <- median(depth_vals)
  z_sd     <- sd(depth_vals)
  z_cv     <- z_sd / z_mean
  
  z_skew   <- e1071::skewness(depth_vals, type = 2)
  z_kurt   <- e1071::kurtosis(depth_vals, type = 2)
  
  # Depth percentiles
  depth_pctls <- quantile(
    depth_vals,
    probs = c(0.10, 0.25, 0.50, 0.75, 0.90),
    names = TRUE
  )
  
  # --- 4. Volume -------------------------------------------------------------
  # Volume = sum of (depth * cell_area) for each cell
  volume_m3 <- sum(depth_vals * cell_area_m2)
  
  # Verify z_mean from volume: should equal volume / area
  z_mean_vol <- volume_m3 / total_area
  
  # --- 5. Morphometric ratios ------------------------------------------------
  # Volume development (Vd): ratio of actual volume to volume of a cone
  # Cone volume = (1/3) * area * z_max
  # Vd = 3 * z_mean / z_max
  # Vd < 1 = concave basin; Vd > 1 = convex; Vd = 1 = conical
  v_d <- 3 * z_mean / z_max
  
  # Relative depth (Zr): max depth relative to lake size (%)
  # Zr = 50 * z_max * sqrt(pi) / sqrt(area)
  z_r <- (50 * z_max * sqrt(pi)) / sqrt(total_area)
  
  # Dynamic ratio: proxy for thermal mixing potential
  # DR = sqrt(area_km2) / z_mean
  area_km2    <- total_area / 1e6
  dynamic_ratio <- sqrt(area_km2) / z_mean
  
  # --- 6. Shoreline & shape metrics ------------------------------------------
  # Derive shoreline from raster if not provided
  if (is.null(shoreline)) {
    lake_mask      <- !is.na(depth_rast)
    shoreline <- terra::as.polygons(lake_mask) |>
      (\(x) x[terra::values(x)[, 1] == 1, ])()   # keep lake pixels only
  }
  
  # Ensure sf for perimeter / area calculations
  shoreline_sf   <- sf::st_as_sf(shoreline)
  perimeter_m    <- sf::st_perimeter(shoreline_sf) |>
    as.numeric() |>
    sum()                           # sum in case of multipolygon
  
  # Shoreline development (DL): ratio of shoreline to circumference of
  # a circle with the same area. DL = 1 for a perfect circle.
  d_l <- perimeter_m / (2 * sqrt(pi * total_area))
  
  # Circularity index (inverse of DL^2 scaled)
  circularity <- (4 * pi * total_area) / perimeter_m^2
  
  # Bounding box elongation ratio (length / width of minimum bounding rectangle)
  bbox        <- sf::st_bbox(shoreline_sf)
  bbox_length <- max(bbox["xmax"] - bbox["xmin"], bbox["ymax"] - bbox["ymin"])
  bbox_width  <- min(bbox["xmax"] - bbox["xmin"], bbox["ymax"] - bbox["ymin"])
  elongation  <- bbox_length / bbox_width
  
  # --- 7. Slope statistics ---------------------------------------------------
  slope_rast  <- terra::terrain(depth_rast, v = "slope", unit = "degrees")
  slope_vals  <- terra::values(slope_rast, na.rm = TRUE) |> as.numeric()
  
  slope_mean   <- mean(slope_vals)
  slope_median <- median(slope_vals)
  slope_sd     <- sd(slope_vals)
  slope_skew   <- e1071::skewness(slope_vals, type = 2)
  
  # Proportion of lake area in slope classes (%)
  pct_slope_gentle   <- mean(slope_vals <  5)  * 100   # < 5°
  pct_slope_moderate <- mean(slope_vals >= 5  & slope_vals < 15) * 100
  pct_slope_steep    <- mean(slope_vals >= 15) * 100   # >= 15°
  
  # --- 8. Hypsographic curve & fitted parameters ----------------------------
  # Build hypsographic curve: cumulative area vs depth
  # Normalise depth to [0, 1] for fitting
  hyps <- .calc_hypsographic(depth_vals, cell_area_m2, n_breaks = 100)
  
  # Fit power model: A(z) = A0 * (1 - z/z_max)^b
  # Rearranged for nls: normalised_area ~ (1 - norm_depth)^b
  hyps_fit  <- .fit_hypsographic_power(hyps, z_max)
  hyps_b    <- hyps_fit$b           # basin shape exponent
  hyps_r2   <- hyps_fit$r_squared   # goodness of fit
  
  # Area under normalised hypsographic curve (trapezoid integration)
  # = 1 for a flat basin filling uniformly; < 1 for concave basins
  hyps_auc <- .trapz(hyps$norm_depth, hyps$norm_area)
  
  # --- 9. Depth zone areas ---------------------------------------------------
  # Percentage of surface area in ecologically relevant depth zones
  pct_area_0_2m   <- mean(depth_vals <= 2)  * 100
  pct_area_0_5m   <- mean(depth_vals <= 5)  * 100
  pct_area_0_10m  <- mean(depth_vals <= 10) * 100
  pct_area_deep   <- mean(depth_vals >  10) * 100
  
  # Volume at depth percentiles (cumulative volume above each depth threshold)
  vol_pctls <- .volume_at_depth_percentiles(depth_vals, cell_area_m2,
                                            probs = c(0.25, 0.50, 0.75))
  
  # --- 10. Compile output ----------------------------------------------------
  results <- list(
    
    # -- Basic geometry
    surface_area_m2     = total_area,
    surface_area_km2    = area_km2,
    volume_m3           = volume_m3,
    perimeter_m         = perimeter_m,
    
    # -- Depth statistics
    z_max               = z_max,
    z_mean              = z_mean,
    z_mean_from_volume  = z_mean_vol,
    z_median            = z_median,
    z_sd                = z_sd,
    z_cv                = z_cv,
    z_skewness          = z_skew,
    z_kurtosis          = z_kurt,
    
    # -- Depth percentiles
    z_p10               = depth_pctls[["10%"]],
    z_p25               = depth_pctls[["25%"]],
    z_p50               = depth_pctls[["50%"]],
    z_p75               = depth_pctls[["75%"]],
    z_p90               = depth_pctls[["90%"]],
    
    # -- Morphometric ratios
    volume_development  = v_d,
    relative_depth      = z_r,
    dynamic_ratio       = dynamic_ratio,
    
    # -- Shape metrics
    shoreline_development = d_l,
    circularity           = circularity,
    elongation_ratio      = elongation,
    
    # -- Slope metrics
    slope_mean            = slope_mean,
    slope_median          = slope_median,
    slope_sd              = slope_sd,
    slope_skewness        = slope_skew,
    pct_slope_gentle      = pct_slope_gentle,
    pct_slope_moderate    = pct_slope_moderate,
    pct_slope_steep       = pct_slope_steep,
    
    # -- Hypsographic curve
    hyps_b_exponent       = hyps_b,
    hyps_r_squared        = hyps_r2,
    hyps_auc_normalised   = hyps_auc,
    
    # -- Depth zone areas (%)
    pct_area_0_2m         = pct_area_0_2m,
    pct_area_0_5m         = pct_area_0_5m,
    pct_area_0_10m        = pct_area_0_10m,
    pct_area_gt10m        = pct_area_deep,
    
    # -- Cumulative volume at depth percentiles (m3)
    vol_at_z25pct         = vol_pctls[["p25"]],
    vol_at_z50pct         = vol_pctls[["p50"]],
    vol_at_z75pct         = vol_pctls[["p75"]]
  )
  
  return(results)
}


#' Apply calc_lake_morphometry to a named list of SpatRasters
#'
#' @description
#' Iterates `calc_lake_morphometry()` over a named list of bathymetric rasters
#' and returns results as a single tidy `data.frame` with one row per lake.
#'
#' Parallel execution is supported via the \pkg{parallel} package (base R).
#' When `parallel = TRUE` the function creates a PSOCK cluster with `n_cores`
#' workers, exports the required functions and packages to each worker, then
#' distributes lakes across workers using `parallel::parLapply()`. The cluster
#' is always stopped on exit — including on error — via `on.exit()`, so no
#' manual cleanup is needed.
#'
#' **Choosing a backend**
#'
#' - `parallel = FALSE` (default): single-threaded `lapply()`. Safe in all
#'   contexts including interactive sessions and RMarkdown.
#' - `parallel = TRUE`: spawns `n_cores` PSOCK workers. Recommended when
#'   processing many lakes (roughly > 20) on a machine with multiple cores.
#'   Note that \pkg{terra} SpatRaster objects are serialised to each worker;
#'   very large rasters may make parallelisation slower than sequential
#'   processing due to inter-process data transfer overhead.
#'
#' @param raster_list A named list of `terra::SpatRaster` objects, one per
#'   lake. If unnamed, lakes are labelled `lake_1`, `lake_2`, etc.
#' @param parallel Logical. If `TRUE`, processing is distributed across
#'   multiple CPU cores using `parallel::makeCluster()` /
#'   `parallel::parLapply()`. Default `FALSE`.
#' @param n_cores Integer. Number of parallel workers to use when
#'   `parallel = TRUE`. Defaults to `parallel::detectCores() - 1L`, leaving
#'   one core free for the main R session. Silently capped at the number of
#'   lakes so no idle workers are created.
#' @param ... Additional arguments passed to `calc_lake_morphometry()`.
#'
#' @return A `data.frame` with one row per successfully processed lake and one
#'   column per morphometric metric, plus a leading `lake_id` column. Lakes
#'   that raise an error are dropped with a `warning()` and do not appear in
#'   the output.
#' @export
#' @examples
#' \dontrun{
#'   files     <- list.files("bathymetry/", pattern = "\\.tif$", full.names = TRUE)
#'   rasters   <- stats::setNames(lapply(files, terra::rast),
#'                                tools::file_path_sans_ext(basename(files)))
#'
#'   # Sequential
#'   results   <- batch_lake_morphometry(rasters)
#'
#'   # Parallel — use 4 cores
#'   results_p <- batch_lake_morphometry(rasters, parallel = TRUE, n_cores = 4)
#' }
batch_lake_morphometry <- function(raster_list,
                                   parallel = FALSE,
                                   n_cores  = parallel::detectCores() - 1L,
                                   ...) {
  
  # --- 0. Name unlabelled lists ----------------------------------------------
  if (is.null(names(raster_list))) {
    names(raster_list) <- paste0("lake_", seq_along(raster_list))
  }
  
  # --- 1. Shared processing closure -----------------------------------------
  # Wraps calc_lake_morphometry in error handling; returns a data.frame or NULL
  .process_one <- \(r) {
    tryCatch(
      calc_lake_morphometry(r, ...) |> as.data.frame(),
      error = \(e) {
        warning("Error processing lake: ", conditionMessage(e))
        NULL
      }
    )
  }
  
  # --- 2. Execute: sequential or parallel -----------------------------------
  if (!parallel) {
    
    results <- lapply(raster_list, .process_one)
    
  } else {
    
    # Cap workers at the number of lakes — no point in idle workers
    n_cores <- min(n_cores, length(raster_list))
    n_cores <- max(n_cores, 1L)   # guard against 0 or negative values
    
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    # Each worker needs the processing function and its dependencies
    parallel::clusterExport(
      cl,
      varlist = c("calc_lake_morphometry",
                  ".calc_hypsographic",
                  ".fit_hypsographic_power",
                  ".trapz",
                  ".volume_at_depth_percentiles"),
      envir   = environment()
    )
    
    parallel::clusterEvalQ(cl, {
      library(terra)
      library(e1071)
      library(sf)
    })
    
    results <- parallel::parLapply(cl, raster_list, .process_one)
  }
  
  # --- 3. Combine results ---------------------------------------------------
  Map(\(df, nm) if (!is.null(df)) cbind(lake_id = nm, df),
      results, names(results)) |>
    (\(x) x[!sapply(x, is.null)])() |>
    do.call(what = rbind)
}


# =============================================================================
# Internal helpers (not exported)
# =============================================================================

#' Build hypsographic curve data frame
#' @keywords internal
#' @noRd
.calc_hypsographic <- function(depth_vals, cell_area_m2, n_breaks = 100) {
  
  z_max      <- max(depth_vals)
  depth_cuts <- seq(0, z_max, length.out = n_breaks + 1)
  
  # Area at or shallower than each depth threshold
  area_at_depth <- vapply(
    depth_cuts,
    \(d) sum(depth_vals <= d) * cell_area_m2,
    FUN.VALUE = numeric(1)
  )
  
  total_area <- max(area_at_depth)
  
  data.frame(
    depth      = depth_cuts,
    area_m2    = area_at_depth,
    norm_depth = depth_cuts / z_max,
    norm_area  = area_at_depth / total_area
  )
}


#' Fit power model to normalised hypsographic curve
#' @keywords internal
#' @noRd
.fit_hypsographic_power <- function(hyps, z_max) {
  
  # Model: norm_area ~ (1 - norm_depth)^b
  # Use nls with a safe starting value
  fit <- tryCatch(
    nls(
      norm_area ~ (1 - norm_depth)^b,
      data    = hyps,
      start   = list(b = 1),
      control = nls.control(maxiter = 200)
    ),
    error = \(e) NULL
  )
  
  if (is.null(fit)) {
    return(list(b = NA_real_, r_squared = NA_real_))
  }
  
  fitted_vals <- fitted(fit)
  ss_res      <- sum((hyps$norm_area - fitted_vals)^2)
  ss_tot      <- sum((hyps$norm_area - mean(hyps$norm_area))^2)
  r2          <- 1 - ss_res / ss_tot
  
  list(
    b         = coef(fit)[["b"]],
    r_squared = r2
  )
}


#' Trapezoidal integration
#' @keywords internal
#' @noRd
.trapz <- function(x, y) {
  dx <- diff(x)
  sum(dx * (utils::head(y, -1) + utils::tail(y, -1)) / 2)
}


#' Cumulative volume at depth percentile thresholds
#' @keywords internal
#' @noRd
.volume_at_depth_percentiles <- function(depth_vals, cell_area_m2,
                                         probs = c(0.25, 0.50, 0.75)) {
  thresholds <- quantile(depth_vals, probs = probs)
  
  vol <- vapply(
    thresholds,
    \(thresh) sum(depth_vals[depth_vals <= thresh]) * cell_area_m2,
    FUN.VALUE = numeric(1)
  )
  
  setNames(as.list(vol), paste0("p", probs * 100))
}

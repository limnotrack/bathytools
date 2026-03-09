# Calculate comprehensive morphometric statistics from a lake bathymetry raster

Computes a comprehensive set of lake morphometric statistics from a
bathymetric raster. Metrics span six categories: basic geometry, depth
statistics, morphometric ratios, shoreline shape, basin slope, and
hypsographic curve parameters. Together these characterise the 3-D form
of the lake basin and are suitable as features for predictive modelling
or comparative limnological analysis.

**1. Basic Geometry**

- `surface_area_m2` / `surface_area_km2`: Planimetric area of the water
  surface, computed as the count of non-NA raster cells multiplied by
  the mean cell area. Units: m² and km².

- `volume_m3`: Total lake volume computed by summing the product of
  depth and cell area across all cells: \$\$V = \sum\_{i=1}^{n} z_i
  \cdot a\$\$ where \\z_i\\ is the depth of cell \\i\\ and \\a\\ is the
  cell area (m²).

- `perimeter_m`: Length of the lake shoreline polygon in metres, derived
  either from the supplied `shoreline_vect` or vectorised from the
  raster extent.

**2. Depth Statistics**

All depth statistics are computed from the vector of per-cell depth
values.

- `z_max`: Maximum depth (m). The deepest recorded point in the basin.

- `z_mean`: Arithmetic mean depth (m): \$\$\bar{z} = \frac{1}{n}
  \sum\_{i=1}^{n} z_i\$\$ Also computed independently as \\V / A\\
  (`z_mean_from_volume`) as a consistency check; both values should be
  nearly identical.

- `z_median`: Median depth (m). More robust than the mean when the depth
  distribution is strongly skewed.

- `z_sd`: Standard deviation of depth (m): \$\$s_z = \sqrt{\frac{1}{n-1}
  \sum\_{i=1}^{n} (z_i - \bar{z})^2}\$\$

- `z_cv`: Coefficient of variation of depth (unitless). Normalises
  variability by the mean, enabling comparison across lakes of different
  sizes: \$\$CV_z = \frac{s_z}{\bar{z}}\$\$

- `z_skewness`: Pearson's moment coefficient of skewness (type 2,
  bias-corrected; via
  [`e1071::skewness`](https://rdrr.io/pkg/e1071/man/skewness.html)).
  Positive values indicate that most of the lake area is shallow with a
  long tail of deeper water (common in glacially scoured basins);
  negative values indicate a flat bottom with steep walls. \$\$g_1 =
  \frac{n}{(n-1)(n-2)} \sum\_{i=1}^{n} \left(\frac{z_i -
  \bar{z}}{s_z}\right)^3\$\$

- `z_kurtosis`: Excess kurtosis of the depth distribution (type 2,
  bias-corrected; via
  [`e1071::kurtosis`](https://rdrr.io/pkg/e1071/man/kurtosis.html)).
  High positive kurtosis indicates a peaked distribution (many cells
  near a single dominant depth); values near zero indicate a flat,
  uniform distribution.

- `z_p10`, `z_p25`, `z_p50`, `z_p75`, `z_p90`: Depth percentiles (m) at
  the 10th, 25th, 50th, 75th, and 90th quantiles of the depth
  distribution. Together these provide a non-parametric summary of the
  hypsographic curve without requiring a functional fit.

**3. Morphometric Ratios**

- `volume_development` (\\V_d\\): Compares the actual basin volume to
  the volume of a cone with the same surface area and maximum depth. A
  cone has \\V_d = 1\\; values \< 1 indicate a concave (bowl-shaped)
  basin where volume is concentrated near the surface; values \> 1
  indicate a convex (flat-bottomed) basin. Introduced by Hutchinson
  (1957): \$\$V_d = \frac{3\\\bar{z}}{z\_{max}}\$\$

- `relative_depth` (\\Z_r\\): Expresses maximum depth as a percentage of
  the mean diameter of the lake surface, allowing depth to be compared
  across lakes of very different sizes (Hutchinson 1957): \$\$Z_r =
  \frac{50\\ z\_{max} \sqrt{\pi}}{\sqrt{A}}\$\$ where \\A\\ is surface
  area in m². High values indicate a deep, small-surface lake; low
  values a shallow, large-surface lake. Strongly linked to
  stratification stability.

- `dynamic_ratio` (\\DR\\): A predictor of thermal stratification and
  wind mixing (Håkanson 1982). High values indicate lakes prone to full
  mixing; low values indicate stable stratification: \$\$DR =
  \frac{\sqrt{A\_{km^2}}}{\bar{z}}\$\$ where \\A\_{km^2}\\ is surface
  area in km² and \\\bar{z}\\ is mean depth in metres. Lakes with \\DR
  \> \approx 0.8\\ are generally considered unstratified.

**4. Shoreline & Shape Metrics**

- `shoreline_development` (\\D_L\\): Ratio of the actual shoreline
  length to the circumference of a circle enclosing the same area. A
  perfect circle gives \\D_L = 1\\; larger values indicate greater
  shoreline irregularity and thus greater littoral habitat complexity
  (Hutchinson 1957): \$\$D_L = \frac{L}{2\sqrt{\pi A}}\$\$ where \\L\\
  is perimeter (m) and \\A\\ is surface area (m²).

- `circularity`: The isoperimetric quotient, an alternative shape index
  bounded on \\0, 1\\ where 1 is a perfect circle. Mathematically the
  reciprocal of \\D_L^2\\: \$\$C = \frac{4\pi A}{L^2}\$\$

- `elongation_ratio`: Ratio of the longest to the shortest axis of the
  axis-aligned bounding box of the shoreline polygon. Values near 1
  indicate a compact, equidimensional lake; large values indicate an
  elongated lake, which affects fetch, internal seiching, and
  circulation patterns. \$\$E = \frac{\max(\Delta x,\\ \Delta
  y)}{\min(\Delta x,\\ \Delta y)}\$\$

**5. Basin Slope Metrics**

Slope is computed cell-by-cell from the bathymetric raster using
`terra::terrain(v = "slope", unit = "degrees")`, which applies a Horn
(1981) finite-difference gradient estimator to the eight-cell
neighbourhood.

- `slope_mean`: Mean basin slope (degrees). Steeper mean slopes are
  associated with younger, tectonically active, or glacially carved
  basins and with greater sediment focusing toward the profundal zone.

- `slope_median`: Median basin slope (degrees). Less sensitive than the
  mean to extreme cliff or wall cells at basin margins.

- `slope_sd`: Standard deviation of slope (degrees). High values
  indicate heterogeneous relief with both very flat and very steep
  zones.

- `slope_skewness`: Skewness of the slope distribution. Positive skew is
  typical (most cells are gently sloping, with a tail of steep wall
  cells); near-zero skew suggests a uniformly graded basin.

- `pct_slope_gentle` / `pct_slope_moderate` / `pct_slope_steep`:
  Percentage of total lake area falling into three ecologically relevant
  slope classes: gentle (\< 5°), moderate (5–15°), and steep (≥ 15°).
  These classes broadly correspond to littoral sediment stability zones
  used in habitat mapping (Håkanson 1977).

**6. Hypsographic Curve Parameters**

The hypsographic (or hypsometric) curve describes the cumulative
planimetric area of the lake at or above each depth, and is the
fundamental 2-D descriptor of basin shape.

- `hyps_b_exponent` (\\b\\): Exponent of the power-law model fitted to
  the normalised hypsographic curve using nonlinear least squares
  ([`stats::nls`](https://rdrr.io/r/stats/nls.html)): \$\$\hat{A}(z) =
  A_0 \left(1 - \frac{z}{z\_{max}}\right)^b\$\$ where \\\hat{A}(z)\\ is
  the area at depth \\z\\ and \\A_0\\ is total surface area. The
  exponent \\b\\ summarises basin concavity in a single parameter: \\b
  \< 1\\ = concave bowl (volume concentrated near surface), \\b = 1\\ =
  conical basin, \\b \> 1\\ = convex or flat- bottomed basin. Equivalent
  to the \\V_d\\ ratio but estimated continuously across all depths.

- `hyps_r_squared`: Coefficient of determination (R²) for the power-law
  fit. Values close to 1 indicate the power model adequately describes
  the basin shape; lower values suggest a more complex or irregular
  hypsograph that may require a higher-order model.

- `hyps_auc_normalised`: Area under the normalised hypsographic curve,
  computed by trapezoidal integration over depth values scaled to \\0,
  1\\. Equivalent to \\\bar{z} / z\_{max}\\ and therefore closely
  related to \\V_d / 3\\. Values near 1 indicate a flat-bottomed basin;
  values near 0.5 indicate a conical basin; values near 0 indicate a
  very concave basin.

**7. Depth Zone Areas**

- `pct_area_0_2m`, `pct_area_0_5m`, `pct_area_0_10m`: Percentage of
  total lake surface area shallower than 2 m, 5 m, and 10 m
  respectively. These thresholds correspond approximately to the
  macrophyte colonisation limit, the euphotic zone in turbid lakes, and
  the epilimnion depth in small temperate lakes.

- `pct_area_gt10m`: Percentage of area deeper than 10 m (profundal zone
  proxy).

- `vol_at_z25pct`, `vol_at_z50pct`, `vol_at_z75pct`: Cumulative water
  volume (m³) contained in the shallowest 25%, 50%, and 75% of depth
  values respectively. These complement the depth percentiles by
  weighting the depth distribution by cell area and reveal how unevenly
  volume is distributed through the water column.

## Usage

``` r
calc_lake_morphometry(
  bathy_raster,
  water_surface_elev = NULL,
  shoreline = NULL,
  depth_positive = FALSE
)
```

## Arguments

- bathy_raster:

  A terra::SpatRaster of lake bathymetry. Depth values should be
  positive (deeper = larger value) or negative (deeper = more negative).
  NA values are treated as outside the lake boundary.

- water_surface_elev:

  Numeric. Elevation of the water surface. Used to compute depth if
  raster contains elevation rather than depth. If NULL, raster values
  are used directly as depths.

- shoreline:

  Optional terra::SpatVector or sf object of the lake shoreline polygon.
  If NULL, the boundary is derived from non-NA raster cells.

- depth_positive:

  Logical. If TRUE (default), deeper water = larger values. Set FALSE if
  your raster uses negative values for depth.

## Value

A named list of morphometric statistics

## References

Håkanson, L. (1977). The influence of wind, fetch, and water depth on
the distribution of sediments in Lake Vanern, Sweden. *Canadian Journal
of Earth Sciences*, 14(3), 397–412.

Håkanson, L. (1981). *A Manual of Lake Morphometry*. Springer-Verlag,
Berlin.

Håkanson, L. (1982). Lake bottom dynamics and morphometry: The dynamic
ratio. *Water Resources Research*, 18(5), 1444–1450.

Horn, B. K. P. (1981). Hill shading and the reflectance map.
*Proceedings of the IEEE*, 69(1), 14–47.

Hutchinson, G. E. (1957). *A Treatise on Limnology. Vol. 1: Geography,
Physics and Chemistry*. Wiley, New York.#'

## Examples

``` r
if (FALSE) { # \dontrun{
  bathy <- terra::rast("my_lake.tif")
  metrics <- calc_lake_morphometry(bathy)
} # }
```

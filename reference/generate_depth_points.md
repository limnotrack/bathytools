# Generate depth points for interpolation

Generate depth points for interpolation

## Usage

``` r
generate_depth_points(
  shoreline,
  islands = NULL,
  depth_points = NULL,
  contours = NULL,
  res = 2,
  subsample = TRUE,
  crs = NULL
)
```

## Arguments

- shoreline:

  sf object of lake shoreline.

- islands:

  sf object of lake islands if present.

- depth_points:

  sf object of depth points or a dataframe of points with columns 'lon'
  and 'lat'. Must contain a "depth" column. If NULL, then contours must
  be provided. If NULL, then contours must be provided.

- contours:

  sf object of depth contours. If NULL, then points must be provided.

- res:

  numeric resolution of output raster in metres.

- subsample:

  logical subsample points to 10000 if TRUE. Default is TRUE.

- crs:

  target coordinate reference system: object of class `crs`, or input
  string for
  [st_crs](https://r-spatial.github.io/sf/reference/st_crs.html)

## Value

sf object of depth points.

## Examples

``` r
shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
package = "bathytools"))
depth_points <- readRDS(system.file("extdata/depth_points.rds",
package = "bathytools"))

# Generate depth points with points ----
depth_points <- generate_depth_points(shoreline = shoreline,
depth_points = depth_points)
#> Generating depth points... [2026-03-09 21:21:01]
#> Finished! [2026-03-09 21:21:01]

# Generate depth points with contours ----
contours <- readRDS(system.file("extdata/depth_contours.rds",
package = "bathytools"))
depth_points <- generate_depth_points(shoreline = shoreline,
contours = contours)
#> Generating depth points... [2026-03-09 21:21:01]
#> Warning: repeating attributes for all sub-geometries for which they may not be constant
#> Warning: repeating attributes for all sub-geometries for which they may not be constant
#> Warning: repeating attributes for all sub-geometries for which they may not be constant
#> Warning: repeating attributes for all sub-geometries for which they may not be constant
#> Warning: repeating attributes for all sub-geometries for which they may not be constant
#> Warning: repeating attributes for all sub-geometries for which they may not be constant
#> Warning: repeating attributes for all sub-geometries for which they may not be constant
#> Warning: repeating attributes for all sub-geometries for which they may not be constant
#> Warning: repeating attributes for all sub-geometries for which they may not be constant
#> Warning: repeating attributes for all sub-geometries for which they may not be constant
#> Warning: repeating attributes for all sub-geometries for which they may not be constant
#> Warning: repeating attributes for all sub-geometries for which they may not be constant
#> Warning: large number of points for interpolation (100695)
#> Finished! [2026-03-09 21:21:05]
```

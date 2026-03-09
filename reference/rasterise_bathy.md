# Generate depth points for interpolation

Generate depth points for interpolation

## Usage

``` r
rasterise_bathy(
  shoreline,
  islands = NULL,
  depth_points = NULL,
  contours = NULL,
  res = 2,
  subsample = TRUE,
  crs,
  method = "nn",
  print_plot = TRUE,
  n = 1,
  m = 1,
  h = 8
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

- method:

  character interpolation method. Options are nearest neighbour ('nn')
  (default), multilevel B-splines ('MBA'), thin plate spline ('tps'),
  and inverse distance weighting ('idw').

- print_plot:

  logical print plot of interpolated raster.

- n:

  initial size of the spline space in the hierarchical construction
  along the x axis. If the rectangular domain is a square, n = m = 1 is
  recommended. If the x axis is k times the length of the y axis, n = 1,
  m = k is recommended. The default is n = 1.

- m:

  initial size of the spline space in the hierarchical construction
  along the y axis. If the y axis is k times the length of the x axis, m
  = 1, n = k is recommended. The default is m = 1.

- h:

  Number of levels in the hierarchical construction. If, e.g., n = m = 1
  and h = 8, the resulting spline surface has a coefficient grid of size
  \\2^h\\ + 3 = 259 in each direction of the spline surface. See
  references for additional information.

## Value

sf object of depth points.

## Examples

``` r
shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
package = "bathytools"))
depth_points <- readRDS(system.file("extdata/depth_points.rds",
package = "bathytools"))
bathy <- rasterise_bathy(shoreline = shoreline, depth_points = depth_points,
crs = 2193)
#> Generating depth points... [2026-03-09 22:35:50]
#> Finished! [2026-03-09 22:35:50]
#> Interpolating to raster... [2026-03-09 22:35:50]
#> Adjusting depths >= 0 to  -0.81 m
#> Finished! [2026-03-09 22:37:34]

```

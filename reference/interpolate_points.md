# Interpolate point depth data to a raster

Interpolate point depth data to a raster

## Usage

``` r
interpolate_points(
  depth_points,
  shoreline,
  islands = NULL,
  crs,
  res = 2,
  method = "nn",
  n = 1,
  m = 1,
  h = 8,
  print_plot = TRUE
)
```

## Arguments

- depth_points:

  sf object of depth points or a dataframe of points with columns 'lon'
  and 'lat'. Must contain a "depth" column. If NULL, then contours must
  be provided. If NULL, then contours must be provided.

- shoreline:

  sf object of lake shoreline.

- islands:

  sf object of lake islands if present.

- crs:

  target coordinate reference system: object of class `crs`, or input
  string for
  [st_crs](https://r-spatial.github.io/sf/reference/st_crs.html)

- res:

  numeric resolution of output raster in metres.

- method:

  character interpolation method. Options are c('MBA', 'tps', 'nn',
  'idw'). Default is 'nn' (nearest neighbor). 'MBA' uses the
  MBA::mba.surf function, 'tps' uses a thin plate spline from the fields
  package, 'nn' uses nearest neighbor interpolation from the terra
  package, and 'idw' uses inverse distance weighting from the terra
  package.

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

- print_plot:

  logical print plot of interpolated raster.

## Value

SpatRaster

## Examples

``` r
shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
package = "bathytools"))
depth_points <- readRDS(system.file("extdata/depth_points.rds",
                                  package = "bathytools"))
depth_points <- generate_depth_points(shoreline = shoreline,
depth_points = depth_points)
#> Generating depth points... [2026-03-09 22:32:16]
#> Finished! [2026-03-09 22:32:16]
bathy_raster <- interpolate_points(depth_points = depth_points, shoreline = shoreline,
crs = 2193)
#> Interpolating to raster... [2026-03-09 22:32:16]
#> Adjusting depths >= 0 to  -0.81 m
#> Finished! [2026-03-09 22:34:01]

```

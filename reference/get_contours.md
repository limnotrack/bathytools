# Get contours from a bathymetric raster

Get contours from a bathymetric raster

## Usage

``` r
get_contours(bathy_raster, surface = 0, depths = NULL)
```

## Arguments

- bathy_raster:

  SpatRaster object with the bathymetry data.

- surface:

  numeric. The surface elevation of the lake. Default is 0.

- depths:

  numeric. The depths at which to calculate the area. If a single
  numeric value is provided, the function will calculate the area at
  each depth from the surface to the minimum depth of the bathymetry
  raster at intervals of the provided value. If a vector of numeric
  values is provided, the function will calculate the area at each depth
  specified in the vector. if `NULL`, the function will use the default
  depths from the model layer structure. Default is `NULL`.

## Value

sf LINESTRING of the contours.

## Examples

``` r
shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
package = "bathytools"))
depth_points <- readRDS(system.file("extdata/depth_points.rds",
package = "bathytools"))
bathy_raster <- rasterise_bathy(shoreline = shoreline,
depth_points = depth_points, crs = 2193, res = 8)
#> Generating depth points... [2026-03-09 21:21:05]
#> Finished! [2026-03-09 21:21:05]
#> Interpolating to raster... [2026-03-09 21:21:05]
#> Adjusting depths >= 0 to  -0.81 m
#> Finished! [2026-03-09 21:21:23]

contours <- get_contours(bathy_raster = bathy_raster)
```

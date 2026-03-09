# Convert bathymetry raster to hypsograph

This function converts a bathymetry raster to a hypsograph, which is a
data.frame with the area of the lake at different depths. The function
takes a bathymetry raster as input and returns a data frame with the
depth and area at each depth.

## Usage

``` r
bathy_to_hypso(bathy_raster, surface = 0, depths = NULL)
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

data.frame

## Examples

``` r
shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
package = "bathytools"))
depth_points <- readRDS(system.file("extdata/depth_points.rds",
package = "bathytools"))
bathy_raster <- rasterise_bathy(shoreline = shoreline,
depth_points = depth_points, crs = 2193)
#> Generating depth points... [2026-03-09 22:28:00]
#> Finished! [2026-03-09 22:28:00]
#> Interpolating to raster... [2026-03-09 22:28:00]
#> Adjusting depths >= 0 to  -0.81 m
#> Finished! [2026-03-09 22:29:47]

hyps <- bathy_to_hypso(bathy_raster = bathy_raster)
```

# Get the shoreline from a bathymetric raster

Get the shoreline from a bathymetric raster

## Usage

``` r
get_shoreline(bathy_raster, dTolerance = NULL)
```

## Arguments

- bathy_raster:

  SpatRaster object with the bathymetry data.

- dTolerance:

  numeric. The tolerance for simplifying the shoreline. Default is the
  maximum resolution of the bathymetric raster. If set to 0, the
  shoreline will not be simplified.

## Value

sf POLYGON of the shoreline.

## Examples

``` r
shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
package = "bathytools"))
depth_points <- readRDS(system.file("extdata/depth_points.rds",
package = "bathytools"))
bathy_raster <- rasterise_bathy(shoreline = shoreline,
depth_points = depth_points, crs = 2193, res = 8)
#> ℹ Generating depth points for interpolation
#> ✔ Generating depth points for interpolation [194ms]
#> 
#> ℹ Interpolating depth points to raster
#> ℹ Adjusting depths >= 0 to -0.81m
#> ℹ Interpolating depth points to raster

#> ✔ Interpolating depth points to raster [18.1s]
#> 
shoreline2 <- get_shoreline(bathy_raster = bathy_raster)
```

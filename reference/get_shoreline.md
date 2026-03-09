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
#> Generating depth points... [2026-03-09 22:32:00]
#> Finished! [2026-03-09 22:32:00]
#> Interpolating to raster... [2026-03-09 22:32:00]
#> Adjusting depths >= 0 to  -0.81 m
#> Finished! [2026-03-09 22:32:15]

shoreline2 <- get_shoreline(bathy_raster = bathy_raster)
```

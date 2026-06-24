# Calculate the volume of a lake using bathymetry data or a hypsograph

Calculate the volume of a lake using bathymetry data or a hypsograph

## Usage

``` r
calculate_lake_volume(
  bathy_raster = NULL,
  hyps = NULL,
  depth = NULL,
  return_rast = FALSE
)
```

## Arguments

- bathy_raster:

  SpatRaster object with the bathymetry data.

- hyps:

  data.frame with columns 'depth' and 'area'.

- depth:

  numeric. The depth to which to calculate the volume. If provided, the
  volume will be calculated to this depth. If not provided, the volume
  will be calculated to the maximum depth of the bathymetry raster or
  the hypsograph.

- return_rast:

  logical. If TRUE, return a raster with the calculated volume in each
  grid cell. Default is FALSE.

## Value

numeric. The volume of the lake in cubic meters (m^3). If `return_rast`
is TRUE, a SpatRaster object with the volume in each grid cell is
returned.

## Examples

``` r
shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
package = "bathytools"))
depth_points <- readRDS(system.file("extdata/depth_points.rds",
package = "bathytools"))
bathy_raster <- rasterise_bathy(shoreline = shoreline,
depth_points = depth_points, crs = 2193)
#> No islands found.
#> ℹ Generating depth points for interpolation
#> ✔ Generating depth points for interpolation [184ms]
#> 
#> ℹ Interpolating depth points to raster
#> ℹ Adjusting depths >= 0 to -0.82m
#> ℹ Interpolating depth points to raster

#> ✔ Interpolating depth points to raster [14.7s]
#> 
calculate_lake_volume(bathy_raster = bathy_raster)
#> [1] 435826452
```

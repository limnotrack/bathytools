# Get lake surface elevation from DEM

Get lake surface elevation from DEM

## Usage

``` r
get_lake_surface_elevation(dem_raster, shoreline)
```

## Arguments

- dem_raster:

  SpatRaster object with the DEM data.

- shoreline:

  sf object of lake shoreline.

## Value

numeric

## Examples

``` r
shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
package = "bathytools"))
dem_raster <- terra::rast(system.file("extdata/dem_32m.tif",
package = "bathytools"))
lake_elev <- get_lake_surface_elevation(dem_raster = dem_raster,
shoreline = shoreline)
#> Lake surface elevation from DEM: 313.3 m
```

# Merge lake bathymetry raster with a DEM raster

Merge lake bathymetry raster with a DEM raster

## Usage

``` r
merge_bathy_dem(
  shoreline,
  bathy_raster,
  dem_raster,
  transform_shoreline = TRUE,
  catchment = NULL,
  crop_dem_to_catchment = TRUE,
  blend_buffer = 10
)
```

## Arguments

- shoreline:

  sf object of lake shoreline.

- bathy_raster:

  SpatRaster object with the bathymetry data.

- dem_raster:

  SpatRaster object with the DEM data.

- transform_shoreline:

  logical. If TRUE, the shoreline will be transformed to match the CRS
  of the DEM raster.

- catchment:

  sf object of the catchment boundary. If provided, the DEM raster will
  be cropped to the catchment boundary.

- crop_dem_to_catchment:

  logical. If TRUE, the DEM raster will be cropped to the catchment
  boundary. Requires catchment to be provided.

## Value

SpatRaster object with the merged bathymetry and DEM data.

## Examples

``` r
shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
package = "bathytools"))
catchment <- readRDS(system.file("extdata/rotoma_catchment.rds",
                                 package = "bathytools"))
depth_points <- readRDS(system.file("extdata/depth_points.rds",
package = "bathytools"))
bathy_raster <- rasterise_bathy(shoreline = shoreline,
depth_points = depth_points, crs = 2193)
#> Generating depth points... [2026-03-09 22:34:03]
#> Finished! [2026-03-09 22:34:03]
#> Interpolating to raster... [2026-03-09 22:34:03]
#> Adjusting depths >= 0 to  -0.81 m
#> Finished! [2026-03-09 22:35:49]

dem_raster <- terra::rast(system.file("extdata/dem_32m.tif",
package = "bathytools"))
dem_bath <- merge_bathy_dem(shoreline = shoreline, bathy_raster = bathy_raster,
dem_raster = dem_raster, catchment = catchment)
#> Resolutions differ. Resampling bathy_raster to DEM resolution.
#> Lake surface elevation from DEM: 313.3 m
```

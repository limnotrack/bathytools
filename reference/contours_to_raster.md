# Convert lake depth contours to a raster

Convert lake depth contours to a raster

## Usage

``` r
contours_to_raster(
  shoreline,
  contours,
  depth_col = "depth",
  res = 10,
  crs = NULL,
  method = c("idw", "near"),
  idw_power = 2
)
```

## Arguments

- contours:

  An sf object containing depth contour lines or polygons, with a
  numeric depth attribute.

- depth_col:

  Character. Name of the column in `contours` containing depth values.

- res:

  Numeric. Output raster resolution in the CRS units (e.g. metres).

- crs:

  An sf/terra-compatible CRS string, EPSG integer, or NULL to inherit
  from `contours`.

- method:

  Character. Interpolation method passed to
  [`terra::interpIDW`](https://rspatial.github.io/terra/reference/interpIDW.html)
  or `"near"` for nearest-contour assignment. Use `"idw"` (default) for
  smooth interpolation.

- idw_power:

  Numeric. IDW power parameter (default 2).

## Value

A
[`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
of interpolated depth values, masked to the lake boundary derived from
the outermost contour.

## Examples

``` r
if (FALSE) { # \dontrun{
  contours <- sf::st_read("lake_contours.gpkg")
  depth_raster <- contours_to_raster(contours, depth_col = "depth_m", res = 5)
  terra::plot(depth_raster)
} # }
```

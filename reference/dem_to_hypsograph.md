# Extend a hypsograph from a DEM

This function calculates a hypsograph from a bathymetry merged with a
DEM, optionally extending it. to include an area above the lake surface.
It can also calculate the area of the lake at different depths, and
optionally mask the DEM to a shoreline.

## Usage

``` r
dem_to_hypsograph(
  shoreline = NULL,
  dem_bath,
  lake_elev = NULL,
  lake_depth = NULL,
  depths = NULL,
  ext_elev = 0
)
```

## Arguments

- shoreline:

  sf object. The shoreline polygon to mask the DEM to.

- dem_bath:

  SpatRaster object. The DEM with bathymetry data.

- lake_elev:

  numeric. The elevation of the lake surface. If NULL, the function will
  calculate the lake surface elevation from the DEM.

- depths:

  numeric. The depths at which to calculate the area. If a single
  numeric value is provided, the function will calculate the area at
  each depth from the surface to the minimum depth of the bathymetry
  raster at intervals of the provided value. If a vector of numeric
  values is provided, the function will calculate the area at each depth
  specified in the vector. if `NULL`, the function will use the default
  depths from the model layer structure. Default is `NULL`.

- ext_elev:

  numeric. The elevation above the lake surface to extend the
  hypsograph. Default is 0, which means no extension.

## Value

data.frame with columns 'elev', 'area', and 'depth'. The 'elev' column
contains the elevation, the 'area' column contains the area of the lake
at that elevation, and the 'depth' column contains the depth of the lake
at that elevation.

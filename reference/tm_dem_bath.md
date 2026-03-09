# Draw a elevation map with bathymetry and surrounding topography

Draw a elevation map with bathymetry and surrounding topography

## Usage

``` r
tm_dem_bath(dem_bath, lake_elev, ext_elev = 0)
```

## Arguments

- dem_bath:

  SpatRaster object with merged DEM and bathymetry data. Can be
  generated using the
  [merge_bathy_dem](https://limnotrack.github.io/bathytools/reference/merge_bathy_dem.md)
  function.

- lake_elev:

  numeric, elevation of the lake surface. Can be extracted using the
  [get_lake_surface_elevation](https://limnotrack.github.io/bathytools/reference/get_lake_surface_elevation.md)
  function.

- ext_elev:

  numeric, additional elevation to extend the lake surface. This allows
  you to generate a flood map. Default is 0.

# Estimate bathymetry from a shoreline polygon using distance transformation.

This function estimates the bathymetry of a lake from a shoreline
polygon using a distance transformation. The function calculates the
distance from the shoreline polygon to the edge of the lake and uses
this distance to estimate the bathymetry. The maximum depth of the lake
is provided as an input parameter. This is the formula used in the
Hollister and Milstead (2010) paper.

## Usage

``` r
estimate_bathymetry(shoreline, hypsograph = NULL, max_depth, res = 10)
```

## Arguments

- shoreline:

  sf object of lake shoreline.

- hypsograph:

  data.frame. A hypsograph data frame with columns 'depth' and 'area'.
  If provided, the function will use the hypsograph to estimate the
  bathymetry. If not provided, the function will use a linear scaling
  based on the maximum depth.

- max_depth:

  numeric. The maximum depth of the lake.

- res:

  numeric resolution of output raster in metres.

## Value

SpatRaster object with the estimated bathymetry.

## References

Hollister, J. W., W.B. Milstead (2010). Using GIS to Estimate Lake
Volume from Limited Data. Lake and Reservoir Management. 26(3)194-199.
[doi:10.1080/07438141.2010.504321](https://doi.org/10.1080/07438141.2010.504321)

## Examples

``` r
shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
package = "bathytools"))
max_depth <- 81.6
bathy <- estimate_bathymetry(shoreline = shoreline, max_depth = max_depth)
terra::plot(bathy)
```

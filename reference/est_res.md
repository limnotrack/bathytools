# Estimate an appropriate raster resolution for a spatial polygon

Scales resolution based on polygon area using a power-law relationship,
so small areas get fine resolution and large areas get coarser
resolution.

## Usage

``` r
est_res(
  x,
  target_cells = 500,
  min_res = NULL,
  max_res = NULL,
  snap_to = c(2, 4, 8, 16, 20, 24, 32, 40, 64)
)
```

## Arguments

- x:

  An sf polygon or multipolygon object

- target_cells:

  Target number of cells across the shorter dimension (default 500)

- min_res:

  Minimum resolution in CRS units (default NULL = no floor)

- max_res:

  Maximum resolution in CRS units (default NULL = no ceiling)

- snap_to:

  Optional numeric vector of "nice" values to snap to (e.g. c(5, 10, 25,
  50, 100))

## Value

A single numeric value: the estimated resolution in CRS units

## Examples

``` r
est_res(nz_coast)                        # auto-scale
#> Error in est_res(nz_coast): could not find function "est_res"
est_res(nz_coast, target_cells = 200)    # coarser
#> Error in est_res(nz_coast, target_cells = 200): could not find function "est_res"
est_res(nz_coast, snap_to = c(10, 25, 50, 100, 250, 500))
#> Error in est_res(nz_coast, snap_to = c(10, 25, 50, 100, 250, 500)): could not find function "est_res"
```

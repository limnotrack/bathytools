# Detect islands within a lake multipolygon

Islands are encoded as inner rings (holes) in the exterior polygon. This
function extracts those holes and returns them as proper polygons.

## Usage

``` r
detect_islands(x)
```

## Arguments

- x:

  An sf object with POLYGON or MULTIPOLYGON geometry

## Value

An sf object containing island polygons, or NULL if none found

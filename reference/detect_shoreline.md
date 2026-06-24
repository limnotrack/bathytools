# Extract the outer shoreline of a lake multipolygon

Returns the exterior boundary of the lake, ignoring any island holes.
The result is a POLYGON (or MULTIPOLYGON) with no holes.

## Usage

``` r
detect_shoreline(x)
```

## Arguments

- x:

  An sf object with POLYGON or MULTIPOLYGON geometry

## Value

An sf object with the outer shoreline polygon(s)

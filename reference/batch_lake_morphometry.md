# Apply calc_lake_morphometry to a named list of SpatRasters

Iterates
[`calc_lake_morphometry()`](https://limnotrack.github.io/bathytools/reference/calc_lake_morphometry.md)
over a named list of bathymetric rasters and returns results as a single
tidy `data.frame` with one row per lake.

Parallel execution is supported via the parallel package (base R). When
`parallel = TRUE` the function creates a PSOCK cluster with `n_cores`
workers, exports the required functions and packages to each worker,
then distributes lakes across workers using
[`parallel::parLapply()`](https://rdrr.io/r/parallel/clusterApply.html).
The cluster is always stopped on exit — including on error — via
[`on.exit()`](https://rdrr.io/r/base/on.exit.html), so no manual cleanup
is needed.

**Choosing a backend**

- `parallel = FALSE` (default): single-threaded
  [`lapply()`](https://rdrr.io/r/base/lapply.html). Safe in all contexts
  including interactive sessions and RMarkdown.

- `parallel = TRUE`: spawns `n_cores` PSOCK workers. Recommended when
  processing many lakes (roughly \> 20) on a machine with multiple
  cores. Note that terra SpatRaster objects are serialised to each
  worker; very large rasters may make parallelisation slower than
  sequential processing due to inter-process data transfer overhead.

## Usage

``` r
batch_lake_morphometry(
  raster_list,
  parallel = FALSE,
  n_cores = parallel::detectCores() - 1L,
  ...
)
```

## Arguments

- raster_list:

  A named list of
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  objects, one per lake. If unnamed, lakes are labelled `lake_1`,
  `lake_2`, etc.

- parallel:

  Logical. If `TRUE`, processing is distributed across multiple CPU
  cores using
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html)
  /
  [`parallel::parLapply()`](https://rdrr.io/r/parallel/clusterApply.html).
  Default `FALSE`.

- n_cores:

  Integer. Number of parallel workers to use when `parallel = TRUE`.
  Defaults to `parallel::detectCores() - 1L`, leaving one core free for
  the main R session. Silently capped at the number of lakes so no idle
  workers are created.

- ...:

  Additional arguments passed to
  [`calc_lake_morphometry()`](https://limnotrack.github.io/bathytools/reference/calc_lake_morphometry.md).

## Value

A `data.frame` with one row per successfully processed lake and one
column per morphometric metric, plus a leading `lake_id` column. Lakes
that raise an error are dropped with a
[`warning()`](https://rdrr.io/r/base/warning.html) and do not appear in
the output.

## Examples

``` r
if (FALSE) { # \dontrun{
  files     <- list.files("bathymetry/", pattern = "\\.tif$", full.names = TRUE)
  rasters   <- stats::setNames(lapply(files, terra::rast),
                               tools::file_path_sans_ext(basename(files)))

  # Sequential
  results   <- batch_lake_morphometry(rasters)

  # Parallel — use 4 cores
  results_p <- batch_lake_morphometry(rasters, parallel = TRUE, n_cores = 4)
} # }
```

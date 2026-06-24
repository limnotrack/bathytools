# Extract Lake Depth at a Given Point

This function extracts the lake depth at a specific point location from
various data sources: bathymetry raster, point depth data, or contour
data. If multiple data sources are provided, the function will use them
in the following priority order: bathymetry raster \> point data \>
contours.

## Usage

``` r
extract_depth_at_point(
  x,
  bathy_raster = NULL,
  depth_points = NULL,
  shoreline = NULL,
  islands = NULL,
  contours = NULL,
  crs = NULL,
  method = c("bilinear", "simple", "nearest", "idw"),
  max_dist = Inf,
  idw_power = 2,
  n_neighbors = 50
)
```

## Arguments

- x:

  numeric or sf POINT. Either a numeric vector of length 2 containing
  the x and y coordinates (c(x, y)), or an sf POINT object. If numeric
  coordinates are provided, the `crs` parameter must also be specified.

- bathy_raster:

  SpatRaster object with the bathymetry data. Optional if `depth_points`
  or `contours` is provided. If provided along with other data sources,
  this will be used preferentially.

- depth_points:

  sf POINT object with depth data. Must contain a 'depth' column.
  Optional if `bathy_raster` or `contours` is provided. If provided
  along with `contours` (but not `bathy_raster`), this will be used
  preferentially.

- shoreline:

  sf POLYGON or MULTIPOLYGON object representing the shoreline. Optional
  but required if `contours` are provided, as the shoreline is used to
  define the lake boundary for interpolation. Ignored if `contours` is
  not provided.

- islands:

  sf POLYGON or MULTIPOLYGON object representing any islands in the
  lake. Default is NULL. Optional but recommended if islands are
  present, as they can affect the interpolation of depth points and
  contours.

- contours:

  sf LINESTRING or MULTILINESTRING object with contour data. Must
  contain a 'depth' column. Optional if `bathy_raster` or `depth_points`
  is provided. This has the lowest priority if multiple data sources are
  provided.

- crs:

  numeric or character. Coordinate reference system (CRS) of the input
  coordinates if `x` is a numeric vector. Can be an EPSG code or proj4
  string. Required if `x` is numeric; ignored if `x` is an sf object.

- method:

  character. Method to use for extraction. Options are:

  - "bilinear": bilinear interpolation from raster (default for raster)

  - "simple": nearest cell value from raster

  - "nearest": nearest neighbor from point data or contours

  - "idw": inverse distance weighting from point data

  Default is "bilinear" for raster data and "nearest" for point/contour
  data.

- max_dist:

  numeric. Maximum distance (in map units) to search for nearest
  neighbor or for IDW interpolation when using point or contour data.
  Default is Inf (no distance limit).

- idw_power:

  numeric. Power parameter for inverse distance weighting. Default is 2.

- n_neighbors:

  integer. Number of nearest neighbors to use for IDW interpolation.
  Default is 50.

## Value

numeric. The depth value at the specified point. Returns NA if no valid
depth can be extracted (e.g., point is outside the raster extent, or no
nearby points/contours are found within max_dist).

## Details

The function performs boundary checks to verify that the point is within
the extent of the data source before attempting extraction. For raster
data, the function will return NA with a warning if the point is outside
the raster extent. For point and contour data, the function will issue a
warning if the point is outside the data extent but will continue
searching, as features may still be found within the specified max_dist
(though NA will be returned if none are found).

The function uses the cli package to provide informative messages about
the extraction method, boundary checks, and results, making it easier to
understand what the function is doing.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load example data
shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
                                 package = "bathytools"))
depth_points <- readRDS(system.file("extdata/depth_points.rds",
                                  package = "bathytools"))

# Create bathymetry raster
bathy <- rasterise_bathy(shoreline = shoreline, depth_points = depth_points,
                         crs = 2193, res = 8)

# Extract depth at a point from raster
depth <- extract_depth_at_point(x = c(2823700, 6404300),
                                bathy_raster = bathy,
                                crs = 2193)

# Extract depth from point data
depth <- extract_depth_at_point(x = c(2823700, 6404300),
                                depth_points = depth_points,
                                crs = 2193,
                                method = "nearest")

# Extract depth using IDW from point data
depth <- extract_depth_at_point(x = c(2823700, 6404300),
                                depth_points = depth_points,
                                crs = 2193,
                                method = "idw",
                                n_neighbors = 5)

# Extract depth from contours
contours <- get_contours(bathy_raster = bathy)
depth <- extract_depth_at_point(x = c(2823700, 6404300),
                                contours = contours,
                                crs = 2193,
                                method = "nearest")
} # }
```

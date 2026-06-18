test_that("extract depth from bathymetry raster works", {
  shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
                                   package = "bathytools"))
  depth_points <- readRDS(system.file("extdata/depth_points.rds",
                                    package = "bathytools"))
  bathy <- rasterise_bathy(shoreline = shoreline, depth_points = depth_points,
                           crs = 2193, res = 8)
  x <- shoreline |> 
    sf::st_point_on_surface() |> 
    sf::st_coordinates() |> 
    c()
  
  # Test extraction with numeric coordinates
  depth <- extract_depth_at_point(x = x,
                                  bathy_raster = bathy,
                                  crs = 2193)
  testthat::expect_type(depth, "double")
  testthat::expect_true(is.numeric(depth))
  testthat::expect_true(depth <= 0)  # Depth should be negative or zero
  
  # Test extraction with bilinear method
  depth_bilinear <- extract_depth_at_point(x = x,
                                           bathy_raster = bathy,
                                           crs = 2193,
                                           method = "bilinear")
  testthat::expect_type(depth_bilinear, "double")
  
  # Test extraction with simple method
  depth_simple <- extract_depth_at_point(x = x,
                                         bathy_raster = bathy,
                                         crs = 2193,
                                         method = "simple")
  testthat::expect_type(depth_simple, "double")
  
  # Test with sf POINT object
  point_sf <- sf::st_as_sf(data.frame(x = x[1], y = x[2]),
                           coords = c("x", "y"), crs = 2193)
  depth_sf <- extract_depth_at_point(x = point_sf,
                                     bathy_raster = bathy)
  testthat::expect_type(depth_sf, "double")
  testthat::expect_equal(depth_sf, depth_bilinear)
})

test_that("extract depth from point data works", {
  depth_points <- readRDS(system.file("extdata/depth_points.rds",
                                    package = "bathytools"))
  
  # Get a coordinate from the point data
  coords <- depth_points[1, ]
  
  # Test nearest neighbor extraction
  depth_nearest <- extract_depth_at_point(x = c(coords[1, 1], coords[1, 2]),
                                          depth_points = depth_points,
                                          crs = sf::st_crs(depth_points),
                                          method = "nearest")
  testthat::expect_type(depth_nearest, "double")
  testthat::expect_equal(depth_nearest, depth_points$depth[1])
  
  # Test with max_dist constraint that should return NA
  depth_na <- extract_depth_at_point(x = c(coords[1] + 10000, coords[2] + 10000),
                                     depth_points = depth_points,
                                     crs = sf::st_crs(depth_points),
                                     method = "nearest",
                                     max_dist = 100)
  testthat::expect_true(is.na(depth_na))
  
  # Test IDW method
  depth_idw <- extract_depth_at_point(x = c(coords[1] + 10, coords[2] + 10),
                                      depth_points = depth_points,
                                      crs = sf::st_crs(depth_points),
                                      method = "idw",
                                      n_neighbors = 5)
  testthat::expect_type(depth_idw, "double")
  testthat::expect_true(!is.na(depth_idw))
  
  # Test IDW with point at exact location (should return exact value without interpolation)
  depth_idw_exact <- extract_depth_at_point(x = c(coords[1], coords[2]),
                                            depth_points = depth_points,
                                            crs = sf::st_crs(depth_points),
                                            method = "idw",
                                            n_neighbors = 5)
  # When point is exactly at a data point, should return that exact depth
  # Use tight tolerance since distance is exactly 0 (no rounding errors expected)
  testthat::expect_equal(depth_idw_exact, depth_points$depth[1],
                        tolerance = 1e-10)
})

test_that("extract depth from contours works", {
  shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
                                   package = "bathytools"))
  depth_points <- readRDS(system.file("extdata/depth_points.rds",
                                    package = "bathytools"))
  bathy <- rasterise_bathy(shoreline = shoreline, depth_points = depth_points,
                           crs = 2193, res = 8)
  contours <- readRDS(system.file("extdata/depth_contours.rds",
                                  package = "bathytools"))
  
  depth_points_sf <- depth_points |> 
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) |> 
    sf::st_transform(2193)
  
  pnt <- shoreline |> 
    sf::st_sample(size = 1)
  x <- pnt |> 
    sf::st_coordinates() |> 
    c()
  
  # Test extraction from contours
  depth_contour <- extract_depth_at_point(x = x,
                                          contours = contours,
                                          crs = 2193,
                                          method = "nearest")
  testthat::expect_type(depth_contour, "double")
  testthat::expect_true(!is.na(depth_contour))
  testthat::expect_true(depth_contour <= 0)  # Depth should be negative or zero
  
  # Test with very small max_dist - should likely return NA since point is not on contour
  depth_na_small <- extract_depth_at_point(x = x,
                                           contours = contours,
                                           crs = 2193,
                                           method = "nearest",
                                           max_dist = 0.01)
  # Very restrictive distance likely results in NA, but verify it's numeric
  testthat::expect_type(depth_na_small, "double")
  testthat::expect_true(is.na(depth_na_small))
  
  # Test with large max_dist - should return a valid depth
  depth_large_dist <- extract_depth_at_point(x = x,
                                             contours = contours,
                                             crs = 2193,
                                             method = "nearest",
                                             max_dist = 10000)
  testthat::expect_type(depth_large_dist, "double")
  testthat::expect_true(!is.na(depth_large_dist))
})

test_that("input validation works", {
  shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
                                   package = "bathytools"))
  depth_points <- readRDS(system.file("extdata/depth_points.rds",
                                    package = "bathytools"))
  depth_points_sf <- depth_points |> 
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) |> 
    sf::st_transform(2193)
  
  bathy <- rasterise_bathy(shoreline = shoreline, depth_points = depth_points,
                           crs = 2193, res = 8)
  
  # Test error when no data source provided
  testthat::expect_error(
    extract_depth_at_point(x = c(2823700, 6404300), crs = 2193),
    "At least one of"
  )
  
  # Test error when numeric coordinates without crs
  testthat::expect_error(
    extract_depth_at_point(x = c(2823700, 6404300), bathy_raster = bathy),
    "crs' must be specified"
  )
  
  # Test error when invalid coordinate length
  testthat::expect_error(
    extract_depth_at_point(x = c(2823700), bathy_raster = bathy, crs = 2193),
    "vector of length 2"
  )
  
  # Test error when invalid method for raster
  testthat::expect_error(
    extract_depth_at_point(x = c(2823700, 6404300), bathy_raster = bathy,
                          crs = 2193, method = "idw"),
    "method must be 'bilinear' or 'simple'"
  )
  
  # Test error when depth_points lacks depth column
  points_no_depth <- depth_points_sf
  points_no_depth$depth <- NULL
  testthat::expect_error(
    extract_depth_at_point(x = c(2823700, 6404300), depth_points = points_no_depth,
                          crs = 2193),
    "must contain a 'depth' column"
  )
  
  # Test warning when multiple points provided
  multi_point <- sf::st_as_sf(
    data.frame(x = c(2823700, 2823800), y = c(6404300, 6404400)),
    coords = c("x", "y"), crs = 2193
  )
  testthat::expect_warning(
    extract_depth_at_point(x = multi_point, bathy_raster = bathy),
    "Multiple points provided"
  )
})

test_that("CRS transformation works correctly", {
  shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
                                   package = "bathytools"))
  depth_points <- readRDS(system.file("extdata/depth_points.rds",
                                    package = "bathytools"))
  bathy <- rasterise_bathy(shoreline = shoreline, depth_points = depth_points,
                           crs = 2193, res = 8)
  
  # Create point in different CRS (WGS84)
  point_wgs84 <- sf::st_as_sf(data.frame(lon = 176.35, lat = -38.03),
                              coords = c("lon", "lat"), crs = 4326)
  
  # Extract should work with automatic CRS transformation
  depth <- extract_depth_at_point(x = point_wgs84, bathy_raster = bathy)
  testthat::expect_type(depth, "double")
})

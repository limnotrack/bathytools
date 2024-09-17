test_that("interpolation of depth points works", {
  shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
                                   package = "bathytools"))
  point_data <- readRDS(system.file("extdata/depth_points.rds",
                                    package = "bathytools"))
  depth_points <- generate_depth_points(shoreline = shoreline,
                                        point_data = point_data)
  bathy <- interpolate_points(point_data = depth_points, shoreline = shoreline,
                              crs = 2193)
  testthat::expect_s4_class(bathy, "SpatRaster")
})

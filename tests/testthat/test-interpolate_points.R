test_that("interpolation of depth points works", {
  shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
                                   package = "bathytools"))
  depth_points <- readRDS(system.file("extdata/depth_points.rds",
                                    package = "bathytools"))
  depth_points <- generate_depth_points(shoreline = shoreline, 
                                        depth_points = depth_points)
  bathy <- interpolate_points(depth_points = depth_points, res = 10,
                              shoreline = shoreline, crs = 2193)
  testthat::expect_s4_class(bathy, "SpatRaster")
})

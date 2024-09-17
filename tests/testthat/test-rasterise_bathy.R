test_that("generate a bathymetry raster file", {
  shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
                                   package = "bathytools"))
  point_data <- readRDS(system.file("extdata/depth_points.rds",
                                    package = "bathytools"))
  bathy <- rasterise_bathy(shoreline = shoreline, point_data = point_data,
                           crs = 2193)
  testthat::expect_s4_class(bathy, "SpatRaster")
})

test_that("generate a bathymetry raster file", {
  shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
                                   package = "bathytools"))
  depth_points <- readRDS(system.file("extdata/depth_points.rds",
                                    package = "bathytools"))
  bathy <- rasterise_bathy(shoreline = shoreline, depth_points = depth_points,
                           crs = 2193, res = 8)
  testthat::expect_s4_class(bathy, "SpatRaster")
  
  hypsograph <- bathy_to_hypso(bathy_raster = bathy)
  
  mm <- terra::minmax(bathy)
  testthat::expect_true(all(mm < 0))
  
  bathy2 <- estimate_bathymetry(shoreline = shoreline, hypsograph = hypsograph,
                                res = 8)
  bathy_err <- calc_bathy_diff(obs = bathy, pred = bathy2)
  median_err <- median(terra::values(bathy_err), na.rm = TRUE)
  testthat::expect_true(median_err < 1)
  
})

test_that("estimate a bathymetry raster file using max depth", {
  shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
                                   package = "bathytools"))
  max_depth <- 81.6
  bathy <- estimate_bathymetry(shoreline = shoreline, max_depth = max_depth)
  testthat::expect_s4_class(bathy, "SpatRaster")
  mm <- terra::minmax(bathy)
  testthat::expect_true(all(mm < 0))

  depth_points <- readRDS(system.file("extdata/depth_points.rds",
                                    package = "bathytools"))
  bathy_raster <- rasterise_bathy(shoreline = shoreline,
                                  depth_points = depth_points, crs = 2193, res = 8)

  bathy_diff <- calc_bathy_diff(obs = bathy_raster, pred = bathy)
  testthat::expect_s4_class(bathy_diff, "SpatRaster")

})

test_that("generate contours from a bathymetry raster file", {
  shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
                                   package = "bathytools"))
  depth_points <- readRDS(system.file("extdata/depth_points.rds",
                                    package = "bathytools"))
  bathy_raster <- rasterise_bathy(shoreline = shoreline, 
                                  depth_points = depth_points,
                                  crs = 2193, res = 8)
  contours <- get_contours(bathy_raster = bathy_raster)
  testthat::expect_s3_class(contours, "sf")

  shoreline2 <- get_shoreline(bathy_raster = bathy_raster)
  testthat::expect_s3_class(shoreline2, "sf")
})

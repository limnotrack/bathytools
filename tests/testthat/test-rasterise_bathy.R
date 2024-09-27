test_that("generate a bathymetry raster file", {
  shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
                                   package = "bathytools"))
  point_data <- readRDS(system.file("extdata/depth_points.rds",
                                    package = "bathytools"))
  bathy <- rasterise_bathy(shoreline = shoreline, point_data = point_data,
                           crs = 2193, res = 8)
  testthat::expect_s4_class(bathy, "SpatRaster")
  mm <- terra::minmax(bathy)
  testthat::expect_true(all(mm < 0))
})

test_that("generate contours from a bathymetry raster file", {
  shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
                                   package = "bathytools"))
  point_data <- readRDS(system.file("extdata/depth_points.rds",
                                    package = "bathytools"))
  bathy_raster <- rasterise_bathy(shoreline = shoreline, point_data = point_data,
                           crs = 2193, res = 8)
  contours <- get_contours(bathy_raster = bathy_raster)
  testthat::expect_s3_class(contours, "sf")

  shoreline2 <- get_shoreline(bathy_raster = bathy_raster)
  testthat::expect_s3_class(shoreline2, "sf")

  # bathy_raster2 <- rasterise_bathy(shoreline = shoreline2, contours = contours,
  #                                  crs = 2193, res = 8)
  # tm_shape(bathy_raster) + tm_raster() + tm_shape(bathy_raster2) + tm_raster()
  # par(mfrow = c(2, 1))
  # bathy_raster3 <- terra::resample(bathy_raster2, bathy_raster)
  # terra::plot(bathy_raster)
  # terra::plot(bathy_raster2)
  #
  # diff <- bathy_raster - bathy_raster3
  # tm_shape(diff) + tm_raster(style = "cont")
  # diff |>
  #   terra::values() |>
  #   mean(na.rm = TRUE)
  # diff |>
  #   terra::values() |>
  #   hist()
  #   # median(na.rm = TRUE)
  #
  # hyps1 <- bathy_to_hypso(bathy_raster = bathy_raster)
  # hyps2 <- bathy_to_hypso(bathy_raster = bathy_raster2)
  # plot(hyps1$area, hyps1$depth, type = "l")
  # lines(hyps2$area, hyps2$depth, col = "red")

})

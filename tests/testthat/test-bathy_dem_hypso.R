test_that("get lake surface elevation from DEM works", {
  shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
                                   package = "bathytools"))
  dem_raster <- terra::rast(system.file("extdata/dem_32m.tif",
                                        package = "bathytools"))
  lake_elev <- get_lake_surface_elevation(dem_raster = dem_raster,
                                          shoreline = shoreline)
  testthat::expect_true(is.numeric(lake_elev))
})

test_that("can merge bathy with DEM", {
  shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
                                   package = "bathytools"))
  catchment <- readRDS(system.file("extdata/rotoma_catchment.rds",
                                   package = "bathytools"))
  point_data <- readRDS(system.file("extdata/depth_points.rds",
                                    package = "bathytools"))
  bathy_raster <- rasterise_bathy(shoreline = shoreline,
                                  point_data = point_data, crs = 2193)
  dem_raster <- terra::rast(system.file("extdata/dem_32m.tif",
                                        package = "bathytools"))
  dem_bath <- merge_bathy_dem(shoreline = shoreline,
                              bathy_raster = bathy_raster, crop_dem_to_catchment = F,
                              dem_raster = dem_raster, catchment = catchment)
  testthat::expect_s4_class(dem_bath, "SpatRaster")

  lake_elev <- get_lake_surface_elevation(dem_raster = dem_raster,
                                          shoreline = shoreline)

  testthat::expect_true(is.numeric(lake_elev))
  testthat::expect_true(lake_elev == 313.3)

  t1 <- tm_dem_bath(dem_bath = dem_bath, lake_elev = lake_elev)
  testthat::expect_s3_class(t1, "tmap")
})

test_that("can generate hypsograph", {
  shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
                                   package = "bathytools"))
  point_data <- readRDS(system.file("extdata/depth_points.rds",
                                    package = "bathytools"))
  bathy_raster <- rasterise_bathy(shoreline = shoreline,
                                  point_data = point_data, crs = 2193)

  hyps <- bathy_to_hypso(bathy_raster = bathy_raster)
  testthat::expect_true(is.data.frame(hyps))

  shore_area <- shoreline |>
    sf::st_area() |>
    units::drop_units()

  area_diff <- max(hyps$area) - shore_area
  frac <- area_diff / shore_area
  # Test less than 1% difference in surface area calculation
  testthat::expect_true(frac < 0.01)

  vol1 <- calculate_lake_volume(bathy_raster = bathy_raster)
  vol2 <- calculate_lake_volume(hyps = hyps)
  vd <- abs(vol1 - vol2)
  frac2 <- vd / vol1
  testthat::expect_true(frac2 < 0.01)
})

test_that("can merge bathy with DEM", {
  shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
                                   package = "bathytools"))
  catchment <- readRDS(system.file("extdata/rotoma_catchment.rds",
                                   package = "bathytools"))
  point_data <- readRDS(system.file("extdata/depth_points.rds",
                                    package = "bathytools"))
  bathy_raster <- rasterise_bathy(shoreline = shoreline,
                                  point_data = point_data, crs = 2193)

  p1 <- plot_raster_3d(x = bathy_raster, shoreline = shoreline, fact = 8)
  testthat::expect_s3_class(p1, "plotly")

  dem_raster <- terra::rast(system.file("extdata/dem_32m.tif",
                                        package = "bathytools"))
  p2 <- plot_raster_3d(x = dem_raster, fact = 8)
  testthat::expect_s3_class(p2, "plotly")

  dem_bath <- merge_bathy_dem(shoreline = shoreline,
                              bathy_raster = bathy_raster,
                              crop_dem_to_catchment = FALSE,
                              dem_raster = dem_raster, catchment = catchment)
  testthat::expect_s4_class(dem_bath, "SpatRaster")
  lake_elev <- get_lake_surface_elevation(dem_raster = dem_raster,
                                          shoreline = shoreline)
  p3 <- plot_raster_3d(x = dem_bath, fact = 0, shoreline = shoreline,
                       lake_elev = lake_elev, split_lake = TRUE,
                       mask_to_lake = TRUE, add_shoreline = TRUE)
  testthat::expect_s3_class(p3, "plotly")

  testthat::expect_true(is.numeric(lake_elev))
  testthat::expect_true(lake_elev == 313.3)

  t1 <- tm_dem_bath(dem_bath = dem_bath, lake_elev = lake_elev)
  testthat::expect_s3_class(t1, "tmap")
})

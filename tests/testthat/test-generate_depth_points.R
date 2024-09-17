test_that("generate depth points with points", {

  shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
                                   package = "bathytools"))
  point_data <- readRDS(system.file("extdata/depth_points.rds",
                                    package = "bathytools"))
  depth_points <- generate_depth_points(shoreline = shoreline,
                                        point_data = point_data)
  testthat::expect_s3_class(depth_points, "sf")
  testthat::expect_true("depth" %in% names(depth_points))
  testthat::expect_true(all(depth_points$depth <= 0))

})

test_that("generate depth points with contours", {
  shoreline <- readRDS(system.file("extdata/rotoma_shoreline.rds",
                                   package = "bathytools"))
  contours <- readRDS(system.file("extdata/depth_contours.rds",
                                    package = "bathytools"))
  depth_points <- generate_depth_points(shoreline = shoreline,
                                        contours = contours)
  testthat::expect_s3_class(depth_points, "sf")
  testthat::expect_true("depth" %in% names(depth_points))
  testthat::expect_true(all(depth_points$depth <= 0))
})

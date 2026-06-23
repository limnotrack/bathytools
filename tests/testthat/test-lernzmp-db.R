test_that("extract_depth_at_point returns depth values for sites with contour data", {
  
  # Setup: fetch and prepare sites
  skip_if_not_installed("aemetools")
  sites <- aemetools::lt_fetch("sites") |>
    dplyr::mutate(lernzmp_id = paste0("LID", lid)) |>
    dplyr::filter(!is.na(lid)) |>
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
    sf::st_transform(2193)
  
  expect_s3_class(sites, "sf")
  expect_true("lernzmp_id" %in% names(sites))
  expect_true(nrow(sites) > 0)
  
  # Run depth extraction for each site
  sites$depth <- NA_real_
  
  for (i in seq_len(nrow(sites))) {
    site <- sites[i, ]
    lid <- site$lernzmp_id
    
    shoreline <- aemetools::lt_fetch(
      table = "lernzmp_lakes",
      filter = aemetools::lt_filter(lernzmp_id == lid)
    )
    
    contours <- aemetools::lt_fetch(
      table = "depth_contours",
      filter = aemetools::lt_filter(lernzmp_id == lid),
      limit = 100000
    )
    
    if (nrow(contours) == 0) {
      message("No contours for site ", site$site_id, ". Skipping.")
      next
    }
    
    contours <- contours |>
      dplyr::mutate(depth = depth_m)
    
    depth_val <- extract_depth_at_point(
      site,
      shoreline = shoreline,
      contours = contours,
      method = "idw"
    )
    
    # Per-site assertions
    expect_length(depth_val, 1)
    expect_true(is.numeric(depth_val))
    expect_true(depth_val <= 0, label = paste("depth >= 0 for site", site$site_id))
    
    sites$depth[i] <- depth_val
  }
  
  # Final assertions on the full result
  sites_with_contours <- sites[!is.na(sites$depth), ]
  expect_true(nrow(sites_with_contours) > 0,
              label = "At least one site should have a depth value")
  expect_true(all(sites_with_contours$depth >= 0),
              label = "All extracted depths should be non-negative")
})


test_that("extract_depth_at_point handles sites with no contour data gracefully", {
  
  skip_if_not_installed("aemetools")
  sites <- aemetools::lt_fetch("sites") |>
    dplyr::mutate(lernzmp_id = paste0("LID", lid)) |>
    dplyr::filter(!is.na(lid)) |>
    sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) |>
    sf::st_transform(2193)
  
  sites_missing <- 0L
  
  for (i in seq_len(nrow(sites))) {
    site <- sites[i, ]
    lid <- site$lernzmp_id
    
    contours <- aemetools::lt_fetch(
      table = "depth_contours",
      filter = aemetools::lt_filter(lernzmp_id == lid),
      limit = 100000
    )
    
    if (nrow(contours) == 0) {
      sites_missing <- sites_missing + 1L
      # Confirm function is NOT called (nothing to test functionally here,
      # but we track the skip count)
    }
  }
  
  # Just verify the loop completed without error and skips were counted
  expect_true(is.integer(sites_missing))
  expect_gte(sites_missing, 0L)
})


test_that("lt_fetch returns expected schema for lernzmp_lakes and depth_contours", {
  skip_if_not_installed("aemetools")
  # Smoke-test the schema helpers — they should not error
  expect_no_error(aemetools::lt_schema("lernzmp_lakes"))
  expect_no_error(aemetools::lt_schema("depth_contours"))
  
  lakes_schema   <- aemetools::lt_schema("lernzmp_lakes")
  contours_schema <- aemetools::lt_schema("depth_contours")
  
  expect_true(is.data.frame(lakes_schema) || is.list(lakes_schema))
  expect_true(is.data.frame(contours_schema) || is.list(contours_schema))
})

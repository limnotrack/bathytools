#' Guess longitude and latitude columns in a data frame.
#' This function looks for common column names that might represent longitude and latitude.
#' It returns a named vector with the guessed longitude and latitude column names.
#' If it cannot find suitable columns, it throws an error with the available column names.
#' @param df A data frame to inspect for longitude and latitude columns.
#' @return A named vector with 'lon' and 'lat' as names and the
#' guessed column names as values.
#' @noRd
#' @importFrom cli cli_abort
guess_lonlat_cols <- function(df) {
  nms <- names(df)
  lon_candidates <- c("lon", "lng", "long", "longitude", "x")
  lat_candidates <- c("lat", "latitude", "y")
  
  lon_col <- nms[match(TRUE, tolower(nms) %in% lon_candidates)]
  lat_col <- nms[match(TRUE, tolower(nms) %in% lat_candidates)]
  
  if (is.na(lon_col) || is.na(lat_col)) {
    cli::cli_abort(c(
      "x" = "Could not guess longitude and latitude columns.",
      "i" = "Found columns: {paste(nms, collapse = ', ')}"
    ))
  }
  
  c(lon = lon_col, lat = lat_col)
}

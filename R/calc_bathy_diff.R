#' Calculate the difference between observed and predicted bathymetry
#'
#' @param obs SpatRaster object with the observed bathymetry.
#' @param pred SpatRaster object with the predicted bathymetry.
#'
#' @return SpatRaster object with the difference between the observed and
#' predicted bathymetry.
#' @export
#'
#' @importFrom terra resample res

calc_bathy_diff <- function(obs, pred) {

  # Observed bathymetry resolution
  obs_res <- terra::res(obs)
  pred_res <- terra::res(pred)

  # If resolution is different, resample the predicted bathymetry
  if (all(obs_res != pred_res)) {
    pred <- terra::resample(pred, obs)
  }
  pred - obs
}

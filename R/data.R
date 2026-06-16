#' Reference data frame for model layer structure.
#'
#' The reference used to structure model layers within the different models.
#' For the gotm_wet model, this is used to estimate the fractions
#' at different depths. Whereas for the glm_aed and dy_cd models, this is used
#' to define the min and max width of the layers.
#'
#' @format ## `model_layer_structure`
#' A data frame with 191 rows and 3 columns:
#' \describe{
#'   \item{zi}{Interface depth (m)}
#'   \item{h}{Layer thickness (m)}
#'   \item{z}{Layer depth (m)}
#' }
#' @source Package development.
"model_layer_structure"

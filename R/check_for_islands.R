#' Check for islands in a set of polygons
#' 
#' @param x sf POLYGON object
#' @return sf POLYGON object with an additional column `is_island` indicating
#' whether the polygon is an island or not.
#' @export
#' 
#' @importFrom dplyr slice mutate filter bind_rows
#' @importFrom sf st_area st_within
#'

check_for_islands <- function(x) {
  
  shore <- x
  if (nrow(x) == 1) {
    return(shore)
  } else {
    # Set largest polygon to lake_shore
    areas <- shore |> 
      dplyr::mutate(area = sf::st_area(geometry))
    shore_idx <- which.max(areas$area)
    shore <- shore |> 
      dplyr::slice(shore_idx) |>
      dplyr::mutate(is_island = FALSE)
    
    # Check if the others are inside the largest polygon
    pot_islands <- x[-shore_idx, ]
    pot_islands <- pot_islands |>
      dplyr::mutate(is_island = rowSums(sf::st_within(geometry, 
                                                      shore$geometry, 
                                                      sparse = FALSE)) > 0)
    
    # 
    islands <- pot_islands |>
      dplyr::filter(is_island)
    if (nrow(islands) == 0) {
      return(shore)
    } else {
      # If there are islands, add them to the shore
      shore <- dplyr::bind_rows(shore, islands)
      return(shore)
    }
  }
}

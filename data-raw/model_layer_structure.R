
model_layer_structure <- data.frame(
  zi = c(seq(0, 5, 0.2),
         seq(5, 20, 0.5),
         seq(20, 50, 1),
         seq(50, 100, 2),
         seq(100, 200, 5),
         seq(200, 500, 10))
) |>
  dplyr::filter(!duplicated(zi)) |>
  dplyr::mutate(h = diff(c(zi, max(zi) + 10)),
                z = zi + (diff(c(zi, NA))) / 2)

usethis::use_data(model_layer_structure, overwrite = TRUE)

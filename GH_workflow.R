# Build README.md
job::job({
  devtools::build_readme()
}, title = "Build README")

attach(loadNamespace("bathytools"), name = "bathytools_all")

library(usethis)
library(devtools)
load_all()
test()
document()
check()

lat <- -36.8898
lon <- 174.46898
user = "145178"
variable = "2m_temperature"
month = 1:5
year = 2022
site = "test"
path = "test"
dir.create(path)
era5_dataset <-  "reanalysis-era5-land"

# ecmwfr::wf_set_key(user = "tadhg.moore6@gmail.com",
#            key = "f1a26c838d5445fcc36b0e4291bc1053",
#            service = "webapi")
# ecmwfr::wf_set_key(user = "145178",
#            key = "8aebf6a4-35d7-48c0-8240-391eafe555bc",
#            service = "cds")
#
#
# download_era5(lat = lat, lon = lon, year = 2020:2021,
#               user = user, path = path)
#
# met <- AEME::convert_era5(lat = lat, lon = lon, year = 2022, site = site, path = path,
#                           format = "AEME")
met_path <- "era5_download"
dir.create(met_path)
met <- aemetools::get_era5_point(lat = lat, lon = lon,
                                 variables = AEME::era5_ref_table$aeme[-1],
                                 years = 2020:2022, format = "aeme",
                                 parallel = TRUE, ncores = 10)

write.csv(met, here::here("inst", "extdata", "lake",  "data", "meteo.csv"),
          row.names = FALSE, quote = FALSE)

library(usethis)
create_package(path = "bathytools")
use_gpl3_license()
use_readme_rmd()
# use_lifecycle_badge("experimental")
use_pkgdown_github_pages()
use_coverage()
use_github_action()
use_github_action("test-coverage")
# use_data()
usethis::use_data_raw()

use_logo("inst/figures/aeme.png")

use_vignette("intro-aeme", title = "Introduction to AEME")



# CRAN packages
use_package("dplyr")
use_package("sf")
use_package("terra")
use_package("units")

# Suggests
use_package("tmap", type = "Suggests")


# GitHub packages
use_dev_package("aemetools", type = "Suggests",
                remote = "github::limnotrack/aemetools")


job::job({
  devtools::check(error_on = "error")
}, title = "devtools - check")

# Code coverage
job::job({
  covr::codecov(
    quiet = FALSE,
    clean = FALSE
  )
  # covr::codecov()
  covr::package_coverage()
}, title = "Code coverage")

# Build pkgdown
job::job({
  pkgdown::build_site_github_pages(new_process = FALSE, install = TRUE)
}, title = "Build pkgdown")


# R CMD CHECK
rcmdcheck::rcmdcheck()


setwd("inst/extdata/lake/")

# Examples ----
tmpdir <- tempdir()
aeme_dir <- system.file("extdata/lake/", package = "AEME")
# Copy files from package into tempdir
file.copy(aeme_dir, tmpdir, recursive = TRUE)
path <- file.path(tmpdir, "lake")
config <- yaml::read_yaml(file.path(path, "aeme.yaml"))
model_controls <- get_model_controls()
inf_factor = c("glm_aed" = 1, "dy_cd" = 1, "gotm_wet" = 1)
outf_factor = c("glm_aed" = 1, "dy_cd" = 1, "gotm_wet" = 1)
model <- c("gotm_wet")
build_aeme(path = path, config = config, model = model,
               model_controls = model_controls, inf_factor = inf_factor, ext_elev = 5,
               use_bgc = TRUE, use_lw = TRUE)
run_aeme(config = config, model = model, verbose = TRUE, path = path)

unlink("45819_wainamu/", recursive = TRUE, force = TRUE)
config <- yaml::read_yaml("aeme.yaml")
model_controls <- read.csv("model_controls.csv")
inf_factor = c("glm_aed" = 1, "dy_cd" = 1, "gotm_wet" = 1)
outf_factor = c("glm_aed" = 1, "dy_cd" = 1, "gotm_wet" = 1)
model <- c("dy_cd")
model <- c("glm_aed")
build_aeme(path = path, config = config, model = model, model_controls = model_controls,
               inf_factor = inf_factor, ext_elev = 5, use_bgc = T, use_lw = T)
run_aeme(config = config, model = model, verbose = TRUE, path = path)

plot(glmtools::get_surface_height("45819_wainamu/glm_aed/output/output.nc"))
glmtools::plot_var(nc_file = "45819_wainamu/glm_aed/output/output.nc",
                   reference = "bottom")
setwd("../../")

library(hexSticker)
library(ggplot2)

p <- ggplot(aes(x = mpg, y = wt), data = mtcars) + geom_point()
p <- p + theme_void() + theme_transparent()
p <- ggplot() + theme_void() + theme_transparent()

sticker(p, package = "", p_color = "black",
        p_family = "mono",
        # p_fontface = "Courier",
        p_size = 30, s_x=1, s_y=.75, s_width=1.3, s_height=1,
        h_fill = "white", h_color = "black",
        filename = "inst/figures/aeme.png")
# configr::write.config(config, file.path = "test2.yaml", write.type = "yaml")

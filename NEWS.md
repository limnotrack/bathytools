# bathytools 0.1.0

## New features

### Core functions

* New `extract_depth_at_point()` function to extract lake depth at a specific point location from various data sources (bathymetry raster, point depth data, or contour data). Includes:
  * Support for multiple extraction methods: bilinear interpolation, nearest neighbor, and inverse distance weighting (IDW)
  * Automatic data source prioritization (raster > points > contours)
  * Boundary checks to verify points are within data extent
  * CLI messaging for informative feedback during extraction
  * Flexible input: accepts numeric coordinates or sf POINT objects
  * Comprehensive unit tests ensuring robustness

* New `calc_lake_morphometry()` and `batch_lake_morphometry()` functions for calculating lake morphometric properties from bathymetric data

* New `calculate_lake_volume()` function to compute lake volume from bathymetry rasters

* New hypsograph estimation capabilities:
  * `estimate_bathymetry()` can now use hypsograph data to estimate bathymetry
  * Enhanced `bathy_to_hypso()` and `dem_to_hypso()` functions for generating hypsographic curves

### Interpolation and rasterization improvements

* `rasterise_bathy()` enhancements:
  * Added multiple interpolation methods with nearest neighbor (nn) as default
  * Improved handling of islands in lake bathymetry
  * Better depth adjustment for non-negative values
  * CLI progress messages for better user feedback
  
* `interpolate_points()` improvements:
  * Added radius calculation for IDW method
  * Enhanced documentation and error handling
  * Improved CLI messaging

* `generate_depth_points()` updates:
  * Added handler for islands
  * Better progress reporting with CLI messages

### Data handling and utilities

* New utility functions in `utils.R`:
  * Helper functions for guessing latitude/longitude columns in dataframes
  * Improved handling of various input formats

* Enhanced CRS (Coordinate Reference System) checking across multiple functions

* Improved shoreline handling and conversion to vector format for terra package

### Visualization

* New `plot_bathy_3d()` and `plot_raster_3d()` functions for 3D visualization of bathymetry

* Updated to tmap v4 compatibility

* `tm_dem_bath()` function for creating thematic maps

### Data extraction

* `get_contours()` - Extract contour lines from bathymetry
* `get_shoreline()` - Extract shoreline from bathymetry
* `get_depths()` - Extract depth values
* `get_lake_depth()` - Get maximum lake depth
* `get_lake_surface_elevation()` - Extract lake surface elevation
* `extract_ext_elev_polygon()` - Extract external elevation polygon

### Data processing

* `merge_bathy_dem()` - Merge bathymetry with digital elevation model (DEM) data with improved extent matching

* `calc_bathy_diff()` - Calculate difference between bathymetric rasters

* Improved duplicate handling and depth corrections throughout

## Package infrastructure

* Added package logo and pkgdown website configuration
  * Logo design with transparent background
  * Favicon and web manifest files
  * Website URL: https://limnotrack.com/bathytools/
  * GitHub repository: https://github.com/limnotrack/bathytools

* Comprehensive documentation using roxygen2:
  * Detailed function documentation with examples
  * Vignettes for key workflows (e.g., merge-bathy-dem)
  * All NAMESPACE entries generated automatically

* Testing infrastructure:
  * testthat framework (>= 3.0.0)
  * Extensive unit tests for core functions
  * Test coverage monitoring via Codecov
  * Continuous integration with R-CMD-check

* Dependencies:
  * Core spatial packages: terra, sf
  * Visualization: plotly, tmap
  * Interpolation: MBA, e1071
  * User interface: cli
  * Data manipulation: dplyr, units
  * Parallelization: parallel

## Documentation improvements

* Updated README with installation instructions and basic examples

* Improved argument naming consistency across functions

* Enhanced function examples with real-world use cases

* Updated to Roxygen2 version 8.0.0

## Bug fixes and refinements

* Fixed filter bug in depth processing
* Corrected depth sign handling (negative values for depths below surface)
* Improved handling of minimum depth values
* Fixed CRS transformation issues
* Better handling of edge cases and boundary conditions
* Removed redundant code and duplicate functions
* Cleaned up unused imports in NAMESPACE

## Contributors

This version includes contributions from:
* GitHub Copilot (copilot-swe-agent[bot]) - Primary development of `extract_depth_at_point()` function and associated tests
* Tadhg Moore (@tadhg-moore) - Package infrastructure, core bathymetry functions, and overall package design
* Chris McBride - Package co-author

## Acknowledgments

This package was developed by [LimnoTrack](http://limnotrack.com/) as part of the Lake Ecosystem Restoration New Zealand Modelling Platform (LERNZmp) project.

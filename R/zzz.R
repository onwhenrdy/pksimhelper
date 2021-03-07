
.onLoad <- function(libname, pkgname) {
  # install beats als unit
  suppressWarnings(units::install_symbolic_unit("beats", warn = TRUE, dimensionless = TRUE))

  invisible()
}

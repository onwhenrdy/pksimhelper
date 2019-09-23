
# evil global variables
.pkg.env <- new.env(parent = emptyenv())

# install of some basic units that are not available out of the box
.install_beasts <- function() units::install_symbolic_unit("beats", warn = T, dimensionless = TRUE)
assign("had_beats", F, .pkg.env)

.onLoad <- function(libname, pkgname) {

  # install beats als dimensionless unit
  tryCatch(.install_beasts(),
           error = function(e) warning("Could not install beasts as a unit"),
           warning = function(w) .pkg.env$had_beats = T)

  invisible()
}


.onUnload <- function(libname, pkgname) {

  # only remove if the unit has not been registered before the loading of the package
  if (!.pkg.env$had_beats)
    units::remove_symbolic_unit("beats")
}
#' on load initialization of reticulate environment

.onLoad <- function(libname, pkgname) {
  devtools::install_github("cellgeni/sceasy")
}
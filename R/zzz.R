#' on load initialization of reticulate environment
sc <- NULL
sr <- NULL
ad <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference
  sc <<- reticulate::import("scanpy", delay_load = TRUE)
  sr <<- reticulate::import("scanorama", delay_load = TRUE)
  ad <<- reticulate::import("annData", delay_load = TRUE)
}
#' on load initialization of reticulate environment
sc <- NULL
sr <- NULL
ad <- NULL
bbknn <- NULL

#' @export
install_python_packages <- function(){
  reticulate::py_install("scipy")
  reticulate::py_install("pandas")
  reticulate::py_install("numpy")
  reticulate::py_install("scanpy")
  system("pip install scanorama")
  system("pip install anndata")
  system("pip install bbknn")
  system("pip install scipy")
}

.onLoad <- function(libname, pkgname) {
  install_python_packages()
  devtools::install_github("cellgeni/sceasy")
  devtools::install_github('satijalab/seurat-wrappers')
}
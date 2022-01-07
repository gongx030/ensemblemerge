#' on load initialization of reticulate environment
sc <- NULL
sr <- NULL
ad <- NULL
bbknn <- NULL

#' @export
install_python_packages <- function(){
  system("pip3 install scipy")
  system("pip3 install pandas")
  system("pip3 install numpy")
  system("pip3 install scanpy")
  system("pip3 install scanorama")
  system("pip3 install anndata")
  system("pip3 install bbknn")
  system("pip3 install scvi")
}

.onLoad <- function(libname, pkgname) {
  install_python_packages()
  devtools::install_github("cellgeni/sceasy")
  install.packages("rsvd")
  devtools::install_github('satijalab/seurat-wrappers')
}
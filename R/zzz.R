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
}

.onLoad <- function(libname, pkgname) {
}
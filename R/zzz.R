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
  system("pip3 install torch==1.8.0+cpu torchvision==0.9.0+cpu torchaudio==0.8.0 -f https://download.pytorch.org/whl/torch_stable.html")
  system("pip3 install -U scikit-learn")
  system("pip3 install SMILE")
  system("pip3 install torch")
}

.onLoad <- function(libname, pkgname) {
  reticulate::py_config()
  # use superassignment to update global reference
  sc <<- reticulate::import("scanpy", delay_load = TRUE)
  sr <<- reticulate::import("scanorama", delay_load = TRUE)
  ad <<- reticulate::import("anndata", delay_load = TRUE, convert = FALSE)
  bbknn <<- reticulate::import("bbknn", delay_load = TRUE)
}
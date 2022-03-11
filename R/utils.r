.check_method <- function(x){
	available_methods <- c(
		"Seurat", 
		"Scanorama", 
		"Harmony", 
		"Liger", 
		"BBKNN", 
		"Uncorrected", 
		"fastMNN", 
		"scVI"
	)
  ### checking valid parameters ###
  if(!all(x%in% available_methods)){
    stop(sprintf("method must be the following: %s", paste(available_methods, collapse = ", ")))
  }
}

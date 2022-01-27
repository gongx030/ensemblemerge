#' set parameters for integration
#'
#' @importFrom methods new
#'
#' @param method character specifying the integration method i.e "Seurat"
#' @return returns a params S4 class depending on input
#' @export
setParams <- function(methods = "Seurat", ...){

	availableMethods = c(
		"Seurat", 
		"Scanorama", 
		"Harmony", 
		"Liger", 
		"BBKNN", 
		"Uncorrected", 
		"fastMNN", 
		"scVI"
	) # list available availableMethods

  ### checking valid parameters ###
  if(!all(methods %in% availableMethods)){
    stop(sprintf("method must be the following: %s", paste(availableMethods, collapse = ", ")))
  }

	params <- new('ParamsList', lapply(methods, function(method){
	  switch(
 			method,
			"Seurat" = new("SeuratParams", ...),
			"Scanorama" = new("ScanoramaParams", ...),
			"Harmony" = new("HarmonyParams", ...),
			"Liger" = new("LigerParams", ...),
			"BBKNN" = new("BBKNNParams", ...),
			"Uncorrected" = new("UncorrectedParams", ...),
			"fastMNN" = new("FastMNNParams", ...),
			"scVI" = new("scVIParams", ...)
		)
	}))

	if (length(params) == 1L){
		return(params[[1L]])
	}else{
		new('EnsembleMergeParams', constituent = params)
	}
}

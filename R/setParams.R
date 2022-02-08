#' Set parameters for integration
#'
#' @param methods The method names (default: 'Seurat')
#' @param ... Additional arguments
#'
#' @importFrom methods new
#'
#' @return If only one method is used, the output is a parameter object for each method; If more than one methods are used, the output is a `EnsembleMergeParams` object
#'
#' @export
#'
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

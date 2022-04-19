#' Run fastMNN merging
#'
#' @param params a FastMNNParams object
#' @param data a Seurat object
#'
#'
#' @return returns a Seurat object with integrated data
#'
run_fastMNN <- function(params, data){

	if (is.null(data@reductions[[params@pca_name]])){
  	data <- RunPCA(
			object = data, 
			npcs = params@npcs, 
			reduction.name = params@pca_name,
			verbose = FALSE
		)
	}
	

}

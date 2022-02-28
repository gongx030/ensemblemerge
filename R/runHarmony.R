#' Run Harmony functions
#'
#' @param params a HarmonyParams object
#' @param data a Seurat object
#'
#' @importFrom Seurat RunPCA
#'
#' @return returns a Seurat object with integrated data
#'
run_Harmony <- function(params, data){

	if (is.null(data@reductions[[params@pca_name]])){
  	data <- RunPCA(
			object = data, 
			npcs = params@npcs, 
			reduction.name = params@pca_name,
			verbose = FALSE
		)
	}

  data <- harmony::RunHarmony(
		object = data, 
		params@batch, 
		theta = params@theta, 
		plot_convergence = FALSE, 
		max.iter.cluster = params@max_iter_cluster,
		assay.use = params@raw_assay,
		reduction.save = params@name,
		verbose = TRUE
	)

	data

}

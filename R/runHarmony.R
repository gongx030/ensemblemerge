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

  data <- RunPCA(
		object = data, 
		npcs = params@npcs, 
		pc.genes = data@var.genes, 
		reduction.name = params@dimreduc_names[["PCA"]]
	)

  data <- harmony::RunHarmony(
		object = data, 
		params@batch, 
		theta = params@theta_harmony, 
		plot_convergence = FALSE, 
		nclust = params@num_clust, 
		max.iter.cluster = params@max_iter_cluster,
		assay.use = 'RNA',
		reduction.save = params@name
	)

	data

}

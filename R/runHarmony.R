#' Run Harmony functions
#'
#' @importFrom Seurat RunPCA
#'
#' @param params a HarmonyMerge object
#' @param data a data object
#' @return a Seurat object
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

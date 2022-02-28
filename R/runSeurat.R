#' Run Seurat functions
#'
#' Adopted from https://satijalab.org/seurat/articles/integration_introduction.html
#'
#' @param params a SeuratParam object
#' @param data a list of Seurat objects or a SeuratList object
#'
#' @importFrom Seurat FindIntegrationAnchors IntegrateData ScaleData RunPCA
#' @importFrom methods is
#'
#' @return returns a Seurat object with integrated data
#'
run_Seurat <- function(params, data){

	# when k.weight is greater than number of cell in any batches, it will result such error:
	# IntegrateData error: Error in idx[i, ] <- res[[i]][[1]]
	# https://github.com/satijalab/seurat/issues/3930
	n <- sapply(data, ncol) # number of cells per cluster
	params@k.weight <- min(n, params@k.weight)
	sprintf('run_Seurat | setting params@k.weight to %d', params@k.weight) %>% message()

  data <- FindIntegrationAnchors(
		object.list = data, 
		dims = 1:params@npcs,
		verbose = FALSE
	)

 	data <- IntegrateData(
		anchorset = data, 
		dims = 1:params@npcs, 
		k.weight = params@k.weight,
		verbose = FALSE
	)

  data <- ScaleData(object = data, verbose = FALSE)

  data <- RunPCA(
		data, 
		npcs = params@npcs, 
		reduction.name = params@name,
		reduction.key = params@reduction_key,
		verbose = FALSE
	)

	data

}


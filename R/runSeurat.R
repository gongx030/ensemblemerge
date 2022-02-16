#' Run Seurat functions
#'
#' @param params a SeuratParam object
#' @param data a Seurat object
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
	params@k.weight <- min(n)
	sprintf('run_Seurat | setting params@k.weight to %d', params@k.weight) %>% message()

  cell_anchors <- FindIntegrationAnchors(object.list = data, dims = 1:params@npcs)
 	data <- IntegrateData(anchorset = cell_anchors, dims = 1:params@npcs, k.weight = params@k.weight)

  if(params@regressUMI && params@scaling) {
    data <- ScaleData(object = data , vars.to.regress = params@vars.to.regress)
  } else if(params@scaling) { # default option
    data <- ScaleData(object = data)
  }

  data <- RunPCA(data, npcs = params@npcs, reduction.name = params@name)

	data

}


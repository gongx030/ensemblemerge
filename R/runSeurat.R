#' Run Seurat functions
#'
#' Adopted from https://satijalab.org/seurat/articles/integration_introduction.html
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

	features <- VariableFeatures(data)

	d <- SplitObject(data, split.by = params@batch)

	# when k.weight is greater than number of cell in any batches, it will result such error:
	# IntegrateData error: Error in idx[i, ] <- res[[i]][[1]]
	# https://github.com/satijalab/seurat/issues/3930
	#
	n <- sapply(d, ncol) # number of cells per cluster
	params@k.weight <- min(n, params@k.weight)
	sprintf('run_Seurat | setting params@k.weight to %d', params@k.weight) %>% message()

  cell_anchors <- FindIntegrationAnchors(
		object.list = d, 
		dims = 1:params@npcs, 
		anchor.features = features
	)

 	d <- IntegrateData(
		anchorset = cell_anchors, 
		dims = 1:params@npcs, 
		k.weight = params@k.weight
	)

  if(params@regressUMI && params@scaling) {
    d <- ScaleData(object = d, vars.to.regress = params@vars.to.regress)
  } else if(params@scaling) { # default option
    d <- ScaleData(object = d)
  }

  d <- RunPCA(
		d, 
		npcs = params@npcs, 
		reduction.name = params@name,
		reduction.key = params@reduction_key
	)

	d@reductions[[params@name]]@assay.used <- data@active.assay
	data[[params@name]] <- d@reductions[[params@name]]
		
	data

}


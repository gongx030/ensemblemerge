#' Run scVI functions
#'
#' @import Seurat
#' @import SingleCellExperiment
#' @import reticulate
#'
#' @param x SummarizedExpirement, SingleCellExperiment, or Seurat object containing single cell counts matrix
#' @return returns a SummarizedExperiment object of the integrated data
#' @export
run_scVI <- function(params, data){

  scvi <- reticulate::import('scvi')
  anndata <- reticulate::import("anndata")

	features <- VariableFeatures(data)
	adata <- anndata$AnnData(X = t(GetAssayData(data, 'counts')[features, ]), obs = data[[params@batch]])

  scvi$settings$seed <- params@seed
  scvi$model$SCVI$setup_anndata(adata, batch_key = params@batch)
	model <- scvi$model$SCVI(adata, n_latent = params@npcs)
	model$train()
	latent <- model$get_latent_representation()
	rownames(latent) <- colnames(data)
	data[[params@name]] <- CreateDimReducObject(
		embeddings = latent,
		key = sprintf('%s_', params@name), 
		assay = DefaultAssay(data)
	)
	data

}

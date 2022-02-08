#' Run scVI functions
#'
#' @param params a scVIParams object
#' @param data a Seurat object
#'
#' @importFrom Seurat VariableFeatures CreateDimReducObject DefaultAssay
#' @importFrom reticulate import
#'
#' @return returns a Seurat object with integrated data
#'
run_scVI <- function(params, data){

  scvi <- import('scvi')
  anndata <- import("anndata")

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

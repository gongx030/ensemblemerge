#' Run BBKNN integration
#'
#' Adopted from https://scanpy-tutorials.readthedocs.io/en/latest/integrating-data-using-ingest.html#BBKNN
#'
#' @param params a scVIParams object
#' @param data a Seurat object
#'
#' @importFrom reticulate import 
#' @importFrom Seurat RunPCA  CreateDimReducObject DefaultAssay as.Graph as.Neighbor
#'
#' @return returns a Seurat object with integrated data
#'
run_BBKNN <- function(params, data){

  bbknn <- import("bbknn")
	sc <- import("scanpy")
	anndata <- import("anndata")

	adata <- anndata$AnnData(
		X = t(GetAssayData(data, 'data')), 
		obs = data[[params@batch]]
	)

 	sc$tl$pca(adata, n_comps = params@npcs)

	sc$external$pp$bbknn(
		adata, 
		batch_key = params@batch
	)
	sc$tl$umap(adata)
	y <- adata$obsm[['X_umap']]
	rownames(y) <- colnames(data)

	data[[params@umap_name]] <- CreateDimReducObject(
		embeddings = y,
		assay = params@raw_assay,
		key = params@umap_key
	)
	data
}

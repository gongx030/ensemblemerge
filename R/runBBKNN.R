#' Run BBKNN integration
#'
#' Adopted from https://github.com/Teichlab/bbknn
#'
#' @param params a scVIParams object
#' @param data a Seurat object
#'
#' @importFrom reticulate import 
#' @importFrom Seurat RunPCA  CreateDimReducObject DefaultAssay
#'
#' @return returns a Seurat object with integrated data
#'
run_BBKNN <- function(params, data){

  bbknn <- import("bbknn")
	sc <- import("scanpy")
	anndata <- import("anndata")

	adata <- anndata$AnnData(
		X = t(GetAssayData(data, 'scale.data')), 
		obs = data[[params@batch]]
	)

 	sc$tl$pca(adata, n_comps = params@npcs)
	bbknn$bbknn(adata, batch_key = params@batch)

	if (params@ridge_regress){
  	if (params@confounder_key == 'leiden'){
 	   	sc$tl$leiden(adata, resolution = 0.4)
    	bbknn$ridge_regression(adata, batch_key= params@batch, confounder_key = params@confounder_key)
		}else{
    	bbknn$ridge_regression(adata, batch_key = params@batch)
		}
   	sc$tl$pca(adata, n_comps = params@npcs)
   	bbknn$bbknn(adata, batch_key = params@batch)
	}

	latent <- adata$obsm[["X_pca"]]
	rownames(latent) <- colnames(data)
	colnames(latent) <- sprintf('%s_%d', params@reduction_key, 1:params@npcs)

	data[[params@name]] <- CreateDimReducObject(
		embeddings = latent,
		key = params@reduction_key,
		assay = DefaultAssay(data)
	)
	data

}

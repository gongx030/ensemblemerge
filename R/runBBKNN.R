#' Run BBKNN integration
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

  data <- RunPCA(
		object = data, 
		npcs = params@npcs + 1L, 
		pc.genes = data@var.genes, 
		reduction.name = params@dimreduc_names[["PCA"]]
	)

	pca <- data[['pca']]@cell.embeddings
	adata <- anndata$AnnData(X = pca, obs = data[[params@batch]])

	# The weird PCA computation part and replacing it with your original values is unfortunately necessary 
	# due to how AnnData innards operate from a reticulate level. 
 	sc$tl$pca(adata)
	adata$obsm$X_pca <- pca

	bbknn$bbknn(adata, params@batch)

	if (params@ridge_regress){
  	if (params@confounder_key == 'leiden'){
 	   	sc$tl$leiden(adata, resolution = 0.4)
    	bbknn$ridge_regression(adata, batch_key= params@batch, confounder_key = params@confounder_key)
    	sc$tl$pca(adata)
    	bbknn$bbknn(adata, batch_key = params@batch)
		}else{
    	bbknn$ridge_regression(adata, batch_key = params@batch)
    	sc$tl$pca(data)
    	bbknn$bbknn(data, batch_key = params@batch)
		}
	}else{
 		bbknn$bbknn(data, batch_key= params@batch)
	}

	# ncol(adata$obsm[["X_pca"]]) = params@npcs - 1
	latent <- adata$obsm[["X_pca"]]
	rownames(latent) <- colnames(data)

	data[[params@name]] <- CreateDimReducObject(
		embeddings = latent,
		key = sprintf('%s_', params@name), 
		assay = DefaultAssay(data)
	)
	data

}

#' Run Scanorama integration
#' 
#' @param params a ScanoramaParams object
#' @param data a Seurat object
#'
#' @importFrom reticulate import
#' @importFrom Seurat VariableFeatures SplitObject GetAssayData CreateDimReducObject DefaultAssay
#'
#' @return returns a Seurat object with integrated data
#'
run_Scanorama <- function(params, data){

  sc <- import("scanpy")
  sr <- import("scanorama")
  anndata <- import("anndata")
  np <- import("numpy")

	features <- VariableFeatures(data)
	object.list <- SplitObject(data, split.by = params@batch)
	object.list <- lapply(object.list, function(x) x[features, ])

	assay_list <- list()
	gene_list <- list()
	for (i in 1:length(object.list)){
		assay_list[[i]] <- as.matrix(t(GetAssayData(object.list[[i]], "data")))
		gene_list[[i]] <- rownames(object.list[[i]])
	}
	integrated.corrected.data <- sr$correct(assay_list, gene_list, return_dimred = TRUE, return_dense = TRUE, dimred = params@npcs)

	intdimred <- do.call(rbind, integrated.corrected.data[[1]])
	rownames(intdimred) <- unlist(sapply(object.list, colnames))

	data[[params@name]] <- CreateDimReducObject(
		embeddings = intdimred,
		key = sprintf('%s_', params@name), 
		assay = DefaultAssay(data)
	)
	data
}

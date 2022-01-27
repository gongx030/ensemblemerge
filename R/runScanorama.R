#' Run Seurat V3 functions
#'
#' @import reticulate
#' @import Seurat
#'
#' @param data SingleCellExperiment object containing single cell counts matrix
#' @param params ScanoramaParams object generated from setParams(method = "BBKNN") function
#' @return returns a SummarizedExperiment object of the integrated data
#' @export
run_Scanorama <- function(params, data){

  sc <- reticulate::import("scanpy")
  sr <- reticulate::import("scanorama")
  anndata <- reticulate::import("anndata")
  np <- reticulate::import("numpy")

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

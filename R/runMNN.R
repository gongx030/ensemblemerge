#' Run fastMNN merging
#'
#' @param params a FastMNNParams object
#' @param data a Seurat object
#'
#' @importFrom Seurat VariableFeatures SplitObject RunPCA as.SingleCellExperiment CreateDimReducObject
#' @importFrom SingleCellExperiment reducedDim
#'
#' @return returns a Seurat object with integrated data
#'
run_fastMNN <- function(params, data){

  data <- RunPCA(data, npcs = params@npcs, verbose = FALSE)
	
	features <- VariableFeatures(data)
	object.list <- SplitObject(data, split.by = params@batch)
	object.list <- lapply(object.list, as.SingleCellExperiment)
	object.list <- lapply(object.list, function(obj) obj[features, ])
	out <- do.call(what = batchelor::fastMNN, args = c(object.list, list(d = params@npcs, k = params@n_neighbors)))

	out <- out[, colnames(data)]

	data[[params@name]] <- CreateDimReducObject(
		embeddings = reducedDim(out, 'corrected'),
		key = sprintf('%s_', params@name), 
		assay = DefaultAssay(data)
	)
	data

}

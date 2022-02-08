#' Run uncorrected merging
#'
#' @param params a UncorrectedMerge object
#' @param data a data object
#' @importFrom Seurat RunPCA
#'
#' @return returns a Seurat object of the integrated data
#'
run_Uncorrected <- function(params, data){

  data <- RunPCA(data, npcs = params@npcs)
	data
}

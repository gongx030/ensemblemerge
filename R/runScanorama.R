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
  ### convert params to python ###
  py$min_genes = params@min_genes
  py$min_cells = params@min_cells
  py$batch_size = params@batch_size
  py$return_dense = params@return_dense
  py$knn = params@knn
  py$svd_solver = params@svd_solver
  py$npcs = params@npcs
  py$batch = params@batch

  ### run BBKNN integration ###
  filepath = system.file("R/runScanorama.py", package = "ensemblemerge")
  data = sceasy::convertFormat(data, from = "sce", to = "anndata")
  py$adata = data
  Sys.sleep(3)
  source_python(filepath)
  integrated = sceasy::convertFormat("temp.h5ad", from = "anndata", to = "seurat")
  integrated = as.SingleCellExperiment(integrated)
  file.remove("temp.h5ad")
  return(integrated)
}
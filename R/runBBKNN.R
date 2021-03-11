#' Run Seurat V3 functions
#'
#' @import reticulate
#' @import Seurat
#'
#' @param data SingleCellExperiment object containing single cell counts matrix
#' @param params BBKNNParams object generated from setParams(method = "BBKNN") function
#' @return returns a SummarizedExperiment object of the integrated data
#' @export
run_BBKNN <- function(params, data){
  ### convert params to python ###
  py$min_genes = params@min_genes
  py$min_cells = params@min_cells
  py$svd_solver = params@svd_solver
  py$scale_factor = params@scale_factor
  py$npcs = params@npcs
  py$nhvg = params@nhvg
  py$batch = params@batch
  py$norm_data = params@norm_data
  py$scaling = params@scaling
  py$regressUMI = params@regressUMI
  py$norm_method = params@norm_method
  py$save_knn = params@save_knn
  py$copy = params@copy
  py$neighbors_within_batch = params@neighbors_within_batch
  py$approx = params@approx
  py$trim = params@trim
  py$n_neighbors = params@n_neighbors

  ### run BBKNN integration ###
  filepath = system.file("R/runBBKNN.py", package = "ensemblemerge")
  data = suppressWarnings(sceasy::convertFormat(data, from = "sce", to = "anndata"))
  py$adata = data
  Sys.sleep(3)
  source_python(filepath)
  integrated = sceasy::convertFormat("temp.h5ad", from = "anndata", to = "seurat")
  integrated = as.SingleCellExperiment(integrated)
  file.remove("temp.h5ad")
  return(integrated)
}
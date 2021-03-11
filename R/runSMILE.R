#' Run Seurat V3 functions
#'
#' @import reticulate
#' @import Seurat
#'
#' @param data SingleCellExperiment object containing single cell counts matrix
#' @param params SMILEParams object generated from setParams(method = "BBKNN") function
#' @return returns a SummarizedExperiment object of the integrated data
#' @export
run_SMILE <- function(params, data){
  filepath = system.file("R/runSMILE.py", package = "ensemblemerge")
  py$adata = suppressWarnings(sceasy::convertFormat(data, from = "sce", to = "anndata"))
  source_python(filepath)
  integrated = sceasy::convertFormat("temp.h5ad", from = "anndata", to = "seurat")
  return(integrated)
}
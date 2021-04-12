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
  ### load up python environment ###
  reticulate::py_config()
  sc <- reticulate::import("scanpy", delay_load = TRUE)
  sr <- reticulate::import("scanorama", delay_load = TRUE)
  ad <- reticulate::import("anndata", delay_load = TRUE, convert = FALSE)
  bbknn <- reticulate::import("bbknn", delay_load = TRUE)
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
  #filepath = system.file("R/runBBKNN.py", package = "ensemblemerge")
  data = suppressWarnings(sceasy::convertFormat(data, from = "sce", to = "anndata"))
  py$adata = data
  Sys.sleep(10)
  py_run_string("import numpy as np
import bbknn
import scanpy as sc

sc.pp.filter_cells(adata, min_genes=300)
sc.pp.filter_genes(adata, min_cells=5)

sc.pp.log1p(adata)
sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver=svd_solver)
sc.pp.neighbors(adata,n_neighbors=int(n_neighbors), n_pcs=int(npcs))

adata.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat

adata_bbknn = bbknn.bbknn(adata, copy=copy, neighbors_within_batch=int(neighbors_within_batch),
                          approx=approx,trim=int(trim), batch_key = batch, n_pcs = int(npcs), use_faiss=False)

sc.tl.pca(adata_bbknn, svd_solver=svd_solver,n_comps=int(npcs))
adata_bbknn.write(filename = 'temp.h5ad')")
  #source_python(filepath)
  integrated = sceasy::convertFormat("temp.h5ad", from = "anndata", to = "seurat")
  integrated = as.SingleCellExperiment(integrated)
  file.remove("temp.h5ad")
  return(integrated)
}

#' Run Seurat V3 functions
#'
#' @import Seurat
#' @import rliger
#' @import SeuratWrappers
#'
#' @param data SingleCellExpirement object containing single cell counts matrix
#' @param params LigerParams object
#'
#' @return returns a SingleCellExperiment object of the integrated data
#' @export
run_Liger <- function(params, data){
  data <- Seurat::as.Seurat(data, counts = "counts", data = NULL)
  liger = seuratToLiger(data, combined.seurat = TRUE, meta.var = params@batch)
  liger <- normalize(liger)
  liger <- selectGenes(liger)
  liger <- scaleNotCenter(liger)

  ### Run integration ###
  liger <- optimizeALS(liger, k = 20)
  liger <- quantile_norm(liger)
  #integrated <- ligerToSeurat(liger)
  return(liger)
}
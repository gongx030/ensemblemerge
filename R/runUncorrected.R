#' Run Seurat V3 scaling and normalizing
#'
#' @import Seurat
#'
#' @param data SummarizedExpirement object containing single cell counts matrix
#' @param params UncorrectedParams object generated from setParams(method="Uncorrected", ...)
#'
#' @return returns a SummarizedExperiment object of the integrated data
#' @export
run_Uncorrected <- function(params, data){
  data = as.Seurat(data, counts = "counts", data = NULL)
  if(params@norm_data){
    data <- NormalizeData(data, normalization.method = params@norm_method, 
                        scale.factor = params@scale_factor)
  }
  if(params@hvg){
    data <- Seurat::FindVariableFeatures(data, 
                                      selection.method = "vst", nfeatures = params@numHVG, 
                                      verbose = FALSE)
  }
  else{
      data@assays$RNA@var.features <- rownames(data)
  }
  if(params@regressUMI && params@scaling) {
    data <- ScaleData(object = data, vars.to.regress = params@vars_to_regress)  # in case of read count data
  } else if(params@scaling) { # default option
    data <- ScaleData(object = data)
  }
  data <- as.SingleCellExperiment(data)
  return(data)
}
#' Merge SummarizedExperiment object
#'
#' This function merges a SummarizedExperiment object using a method specified
#' by the method parameter and returns the merged object as a SummarizedExperiment
#' object
#'
#' @param x SummarizedExpirement object containing single cell counts matrix
#' @return returns a SummarizedExperiment object of the integrated data
#' @export
setMethod("Merge", "SummarizedExperiment", function(x, split.by = NULL, batch_label = NULL, merge.which = NULL, method="Seurat") {

    x = as(x, "SingleCellExperiment") # set x as SingleCellExperiment

    availableMethods = c("Seurat", "Scanorama", "Harmony", "LIGER") # list available availableMethods

    ### checking valid parameters ###
    if(!(method %in% availableMethods)){
      stop(sprintf("method must be the following: %s", paste(availableMethods, collapse = ", ")))
    }
    if(!(batch_label %in% names(colData(x)))){
      stop(sprintf("must merge by: %s", paste(names(colData(x)), collapse = ", ")))
    }

    ### running integration ###
    x = switch(method,
              "Seurat" = run_Seurat(x, batch_label),
              "Scanorama" = run_Scanorama(x),
              "Harmony" = run_Harmony(x),
              "LIGER" = run_LIGER(x))
    return(x)
})
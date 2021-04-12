#' Run kBET scoring
#'
#' @import SummarizedExperiment
#'
#' @param x SummarizedExperiment object containing an integrated single cell counts matrix
#' @return returns a kBET score object
#' @export
run_kBET <- function(x, batch_label){
  if (!requireNamespace("kBET", quietly = TRUE)) {
    stop("Package kBET needed for this function to work. Please install it.",
      call. = FALSE)
  }
  score = kBET::kBET(t(as.matrix(assays(x)[["logcounts"]])), colData(x)[,"batch"], plot = TRUE)
  return(score)
}
#' Run kBET scoring
#'
#' @import kBET
#'
#' @param x SummarizedExpirement object containing an integrated single cell counts matrix
#' @return returns a kBET score object
#' @export
run_kBET <- function(x, batch_label){
  score = kBET(t(as.matrix(assays(x)[["logcounts"]])), colData(x)[,"batch"], plot = TRUE)
  return(score)
}
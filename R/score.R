#' score SingleCellExperiment object
#'
#' @import SingleCellExperiment
#'
#' @param x SingleCellExperiment object containing an integrated single cell counts matrix
#' @return returns an output score statistic
#' @export
setMethod("Score", "SingleCellExperiment", function(x, batch_label = NULL, method="kBET") {
  score = switch(method, 
                "kBET" = run_kBET(x, batch_label = batch_label),
                "ARI" = run_ARI(x, batch_label = batch_label),
                "LSI" = run_LSI(x, batch_label = batch_label),
                "ASW" = run_ASW(x, batch_label = batch_label))
  return(score)
})
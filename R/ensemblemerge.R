#' The core ensemble function
#'
#' @param params a EnsembleMergeParams object
#' @param data a Seurat object
#' @importFrom methods as
#' @return a list of two elememts: (1) the intergrated cell-cell neighboring graph, and (2) weight of each constituent method
#'
ensemblemerge_core <- function(params, data){

	cn <- colnames(data[[1L]])
  data <- lapply(1:length(params@constituent), function(i) getNeighborGraph(params@constituent[[i]], data[[i]]))
	ng <- lapply(1:length(params@constituent), function(i) data[[i]][[params@constituent[[i]]@snn_name]])
	ng <- lapply(1:length(params@constituent), function(i) ng[[i]][cn, cn])
	agreement <- Reduce('+', lapply(ng, function(x) x > 0)) 
	agreement <- as(agreement, 'dgCMatrix')
	agreement <- agreement / length(params@constituent)

	# kernel target alignment
	wt <- sapply(1:length(ng), function(i) sum(ng[[i]] * agreement) / sqrt(sum(ng[[i]] * ng[[i]]) * sum(agreement * agreement)))	
	wt <- (wt -  min(wt)) / (max(wt) - min(wt))
	sigmoid <- function(x, a = 1, b = 1){1/(1+a*(x/(1-x))^-(b*2))}
	wt <- sigmoid(wt, a = params@sigma_a, b = params@sigma_b)
	wt <- (wt -  min(wt)) / (max(wt) - min(wt))

	ng <- Reduce('+', lapply(1:length(ng), function(i) ng[[i]] * wt[i]))
	ng <- ng / sum(wt)
	list(ng = ng, weight = wt)
}

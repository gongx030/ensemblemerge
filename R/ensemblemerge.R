
ensemblemerge_core <- function(params, data, ...){

	cn <- colnames(data[[1L]])
  data <- lapply(1:length(params@constituent), function(i) getNeighborGraph(params@constituent[[i]], data[[i]]))
	ng <- lapply(1:length(params@constituent), function(i) data[[i]][[params@constituent[[i]]@snn_name]])
	ng <- lapply(1:length(params@constituent), function(i) ng[[i]][cn, cn])
	agreement <- Reduce('+', lapply(ng, function(x) x > 0)) > 1L	# find the agreement edges
	agreement <- as(agreement, 'dgCMatrix')

	# kernel target alignment
	wt <- sapply(1:length(ng), function(i) sum(ng[[i]] * agreement) / sqrt(sum(ng[[i]] * ng[[i]]) * sum(agreement * agreement)))	
	ng <- Reduce('+', lapply(1:length(ng), function(i) ng[[i]] * wt[i]))
	ng <- ng / sum(wt)
	ng
}

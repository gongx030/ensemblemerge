#' @importFrom methods callNextMethod 
#'
setMethod('initialize', 'BaseCluster', function(.Object, check_dependencies = TRUE, ...){
	if (check_dependencies)
		.check_dependences(.Object)
	.Object@cluster_name <- .Object@name
	.Object@snn_name <- sprintf('%sSNN', .Object@name)
	.Object@knn_name <- sprintf('%sKNN', .Object@name)
	callNextMethod(.Object, check_dependencies = check_dependencies, ...)
})


#' The LouvainCluster class
#'
#' @slot k_param Defines k for the k-nearest neighbor algorithm
#'
setClass(
	'LouvainCluster',
	representation(
		snn_name = 'character',
		knn_name = 'character',
		k_param = 'integer'
	),
	contains = c('BaseCluster'),
	prototype(
		name = 'LouvainCluster',
		k_param = 20L	
	)
)

#' Cluster a Seurat object by Louvain method
#'
#' @param x a Seurat object
#' @param params a LouvainCluster object
#' @param ... Additional arguments
#' @return returns a data object with PCA embedding
#' @importFrom methods is
#' @importFrom Seurat FindClusters
#'
setMethod(
	'Cluster',
	signature(
		x = 'Seurat',
		params = 'LouvainCluster'
	),
	function(
		x,
		params,
		...
	){

		# to be implemented
		# 1. whether the embedding is available
		# stopifnot(valid(x, params))	

		x <- x %>%
			FindNeighbors(
				compute.SNN = TRUE,
				reduction = params@embedding@reduction_name,
				dims = 1:params@embedding@ndims,
				k.param = params@k_param,
				graph.name = c(params@knn_name, params@snn_name),
				verbose = FALSE
			)

		x <- x %>% 
			FindClusters(
				graph.name = params@snn_name,
				modularity.fxn = 1,
				resolution = 0.8,
				algorithm = 1,
				random.seed = params@seed,
				verbose = FALSE
			)

		 x[[params@cluster_name]] <- x[['seurat_clusters']]	
		 x[['seurat_clusters']] <- NULL
		 x

	}
)

#' Cluster a SeuratList
#'
#' @param x a SeuratList object
#' @param params a BaseCluster object
#' @param ... Additional arguments
#' @return returns a SeuratList object with clusters
#'
setMethod(
	'Cluster',
	signature(
		x = 'SeuratList',
		params = 'BaseCluster'
	),
	function(
		x,
		params,
		...
	){

		for (i in 1:length(x)){
			x[[i]] <- Cluster(x[[i]], params, ...)
		}
		x
	}
)

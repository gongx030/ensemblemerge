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
				algorithm = 1,	# Louvain
				random.seed = params@seed,
				verbose = FALSE
			)

		 x[[params@cluster_name]] <- x[['seurat_clusters']]	
		 x[['seurat_clusters']] <- NULL
		 x

	}
)


#' The LeidenCluster class
#'
#' @slot k_param Defines k for the k-nearest neighbor algorithm
#'
setClass(
	'LeidenCluster',
	representation(
		snn_name = 'character',
		knn_name = 'character',
		k_param = 'integer'
	),
	contains = c('BaseCluster'),
	prototype(
		name = 'LeidenCluster',
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0'),
			new('PythonPackage', package_name = 'leidenalg', package_version = '0.8.8')
		),
		k_param = 20L	
	)
)

#' Cluster a Seurat object by Leiden method
#'
#' @param x a Seurat object
#' @param params a LeidenCluster object
#' @param ... Additional arguments
#' @return returns a data object with PCA embedding
#' @importFrom methods is
#' @importFrom Seurat FindClusters
#'
setMethod(
	'Cluster',
	signature(
		x = 'Seurat',
		params = 'LeidenCluster'
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
				algorithm = 4, # Leiden
				random.seed = params@seed,
				verbose = FALSE
			)

		 x[[params@cluster_name]] <- x[['seurat_clusters']]	
		 x[['seurat_clusters']] <- NULL
		 x

	}
)

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
		embedding = 'BaseEmbed',
		k_param = 'integer'
	),
	contains = c('BaseCluster'),
	prototype(
		name = 'LouvainCluster',
		k_param = 20L	
	)
)


#' @importFrom methods callNextMethod
#'
setMethod('initialize', 'LouvainCluster', function(.Object, check_dependencies = TRUE, ...){
	.Object <- callNextMethod(.Object, check_dependencies = check_dependencies, ...)
	.Object@snn_name <- sprintf('%sSNN', .Object@name)
	.Object@knn_name <- sprintf('%sKNN', .Object@name)
	callNextMethod(.Object, check_dependencies = check_dependencies, ...)
})

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
		embedding = 'BaseEmbed',
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

#' @importFrom methods callNextMethod
#'
setMethod('initialize', 'LeidenCluster', function(.Object, check_dependencies = TRUE, ...){
	.Object <- callNextMethod(.Object, check_dependencies = check_dependencies, ...)
	.Object@snn_name <- sprintf('%sSNN', .Object@name)
	.Object@knn_name <- sprintf('%sKNN', .Object@name)
	callNextMethod(.Object, check_dependencies = check_dependencies, ...)
})

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


#' The scLCACluster class
#'
setClass(
	'scLCACluster',
	representation(
		clust_max = 'integer',
		training_set_size = 'integer',
		preprocess = 'BasePreprocess',
		zerocorrection = 'numeric',
		cor_thresh = 'numeric'
	),
	contains = c('BaseCluster'),
	prototype(
		name = 'scLCACluster',
		clust_max = 10L,
		training_set_size = 1000L,
		zerocorrection = 0.25,
		cor_thresh = 0.5,
		dependences = list(
			new('RPackage', package_name = 'RMTstat', package_version = '0.3'),
			new('RPackage', package_name = 'mclust', package_version = '5.4.9'),
			new('RPackage', package_name = 'cluster', package_version = '2.1.2'),
			new('RPackage', package_name = 'scLCA', package_version = '0.0.0.9000')
		)
	)
)

#' Cluster a Seurat object by scLCA (https://bitbucket.org/scLCA/single_cell_lca/src/master/)
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
		params = 'scLCACluster'
	),
	function(
		x,
		params,
		...
	){

		# to be implemented
		# 1. whether the embedding is available
		# stopifnot(valid(x, params))	

		raw_assay <- params@preprocess@raw_assay

		batch <- NULL
		if (!is.null(x[[params@preprocess@batch]])){
			if (length(unique(x[[params@preprocess@batch]][, 1])) > 1)
				batch <- x[[params@preprocess@batch]][, 1]
		}

    results <- myscLCA(
			x@assays[[raw_assay]]@data,
			cor.thresh = params@cor_thresh,
			clust.max = params@clust_max, 
			trainingSetSize = params@training_set_size, 
			datBatch = batch,
			outlier.filter = FALSE, 
			zerocorrection = params@zerocorrection
		)

	 x[[params@cluster_name]] <- results[[1]]
	 x

	}
)

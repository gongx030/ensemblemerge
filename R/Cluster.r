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
#' @return returns a data object with clustering results in meta data
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
#' @return returns a data object with clustering results in meta data
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
#' @return returns a data object with clustering results in meta data
#' @importFrom methods is
#' @importFrom Seurat FindClusters
#' @references Cheng C, Easton J, Rosencrance C, Li Y, Ju B, Williams J, Mulder HL, Pang Y, Chen W, Chen X. Latent cellular analysis robustly reveals subtle diversity in large-scale single-cell RNA-seq data. Nucleic Acids Res. 2019 Dec 16;47(22):e143. doi: 10.1093/nar/gkz826. PMID: 31566233; PMCID: PMC6902034.
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

    results <- scLCA::myscLCA(
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


#' The scCCESSKmeansCluster class
#'
setClass(
	'scCCESSCluster',
	representation(
		preprocess = 'BasePreprocess',
		criteria_method = 'character',
		clust_min = 'integer',
		clust_max = 'integer',
		ensemble_sizes = 'integer',
		cores = 'integer',
		batch_size = 'integer',
		max_random_projection = 'integer',
		ndims = 'integer',
		hidden_dims = 'integer',
		learning_rate = 'numeric',
		epochs = 'integer'
	),
	contains = c('BaseCluster'),
	prototype(
		criteria_method = 'NMI',
		clust_min = 5L,
		clust_max = 15L,
		ensemble_sizes = 10L,
		cores = 8L,
		batch_size = 64L,
		max_random_projection = 2048L,
		ndims = 16L,
		hidden_dims = 128L,
		learning_rate = 0.001,
		epochs = 100L,
		dependences = list(
			new('RPackage', package_name = 'tensorflow', package_version = '2.8.0'),
			new('RPackage', package_name = 'keras', package_version = '2.8.0'),
			new('RPackage', package_name = 'clue', package_version = '0.3.60'),
			new('RPackage', package_name = 'scCCESS', package_version = '0.2.0'),
			new('PythonPackage', package_name = 'tensorflow', package_version = '2.8.0'),
			new('PythonPackage', package_name = 'keras', package_version = '2.8.0')
		)
	)
)

setClass(
	'scCCESSKmeansCluster',
	representation(
	),
	contains = 'scCCESSCluster',
	prototype(
		name = 'scCCESSKmeansCluster'
	)	
)

#' Cluster a Seurat object by scCCESS K-means (https://github.com/PYangLab/scCCESS)
#'
#' @param x a Seurat object
#' @param params a LeidenCluster object
#' @param ... Additional arguments
#' @return returns a data object with clustering results in meta data
#' @importFrom methods is
#' @importFrom Seurat FindClusters
#' @importFrom stats kmeans
#' @references Geddes TA, Kim T, Nan L, Burchfield JG, Yang JYH, Tao D, et al. Autoencoder-based cluster ensembles for single-cell RNA-seq data analysis. BMC Bioinformatics. 2019;20:660.
#' According to https://github.com/PYangLab/scCCESS#installation, we used an unofficial `SIMLR` version (https://github.com/yulijia/SIMLR)
#'
setMethod(
	'Cluster',
	signature(
		x = 'Seurat',
		params = 'scCCESSKmeansCluster'
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

		# by default genes_as_rows = T and scale =F 
		# https://github.com/PYangLab/scCCESS/blob/5a79afe0e608717a1a929124b7e7ccf2c66df176/R/scCCESS.R#L330

		k <- scCCESS::estimate_k(
			x@assays[[raw_assay]]@scale.data,
			seed = params@seed, 
			cluster_func = function(x,centers) { 
				set.seed(42)
				kmeans(x, centers)
			},
			criteria_method = params@criteria_method,
			krange = params@clust_min:params@clust_max,
			ensemble_sizes = params@ensemble_sizes,
			cores = params@cores
		)

		cluster <- scCCESS::ensemble_cluster(
			x@assays[[raw_assay]]@scale.data,
			seed = params@seed,
			cluster_func = function(x) {
				set.seed(1)
				kmeans(x, centers = k)
			}, 
			cores = params@cores,
			genes_as_rows = TRUE, 
			ensemble_sizes = params@ensemble_sizes,
			verbose = 0, 	# silent
			scale = FALSE, 
			batch_size = params@batch_size,
			max_random_projection = params@max_random_projection,
			encoded_dim = params@ndims,
			hidden_dims = params@hidden_dims,
			learning_rate = params@learning_rate,
			epochs = params@epochs
		)

		# ensemble_cluster returns a list of length len(ensemble_sizes) containing vectors of
		# consensus clusters per cell. Each ensemble clustering is generated
		# using a number of individual clusterings given by the
		# corresponding element in the ensemble_sizes argument.

	 	x[[params@cluster_name]] <- cluster[[1]]
	 	x

	}
)


setClass(
	'scCCESSSIMLRCluster',
	representation(
	),
	contains = 'scCCESSCluster',
	prototype(
		name = 'scCCESSSIMLRCluster',
		dependences = list(
			new('RPackage', package_name = 'SIMLR', package_version = '1.18.0')
		)
	)	
)

#' Cluster a Seurat object by scCCESS SMILR (https://github.com/PYangLab/scCCESS)
#'
#' @param x a Seurat object
#' @param params a LeidenCluster object
#' @param ... Additional arguments
#' @return returns a data object with clustering results in meta data
#' @importFrom methods is
#' @importFrom Seurat FindClusters
#' @references Geddes TA, Kim T, Nan L, Burchfield JG, Yang JYH, Tao D, et al. Autoencoder-based cluster ensembles for single-cell RNA-seq data analysis. BMC Bioinformatics. 2019;20:660.
#'
setMethod(
	'Cluster',
	signature(
		x = 'Seurat',
		params = 'scCCESSSIMLRCluster'
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

		# by default genes_as_rows = T and scale =F 
		# https://github.com/PYangLab/scCCESS/blob/5a79afe0e608717a1a929124b7e7ccf2c66df176/R/scCCESS.R#L330
		# to be consistent, we used scaled data for 

		k <- scCCESS::estimate_k(
			x@assays[[raw_assay]]@scale.data,
			seed = params@seed, 
			cluster_func = function(x, centers) {
				set.seed(42);
				SIMLR::SIMLR_Large_Scale(t(x), c = centers,kk = 15)
			},
			criteria_method = params@criteria_method,
			krange = params@clust_min:params@clust_max,
			ensemble_sizes = params@ensemble_sizes,
			cores = params@cores
		)

		cluster <- scCCESS::ensemble_cluster(
			x@assays[[raw_assay]]@scale.data,
			seed = params@seed,
			cluster_func = function(x) {
				set.seed(1)
				SIMLR::SIMLR_Large_Scale(t(x), c=k,kk=15)
			}, 
			cores = params@cores,
			genes_as_rows = TRUE, 
			ensemble_sizes = params@ensemble_sizes,
			verbose = 0, 	# silent
			scale = FALSE, 
			batch_size = params@batch_size,
			max_random_projection = params@max_random_projection,
			encoded_dim = params@ndims,
			hidden_dims = params@hidden_dims,
			learning_rate = params@learning_rate,
			epochs = params@epochs
		)

		# ensemble_cluster returns a list of length len(ensemble_sizes) containing vectors of
		# consensus clusters per cell. Each ensemble clustering is generated
		# using a number of individual clusterings given by the
		# corresponding element in the ensemble_sizes argument.

	 	x[[params@cluster_name]] <- cluster[[1]]
	 	x

	}
)

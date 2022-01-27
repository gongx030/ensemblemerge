# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Validity
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @importFrom reticulate py_config py_run_string
#'
check_package <- function(object){

	if (length(object@package_name) == 0){
		return(NULL)

	}else{

		if (length(object@package_version ) == 0){
			txt <- 'package version is empty'
			return(txt)
		}

		if (inherits(object, 'RPackage')){

			is_available <- require(object@package_name, character.only = TRUE)
	
			if (!is_available){
				txt <- sprintf('R package %s is not available', object@package_name)
				return(txt)
			}else{
				if (packageVersion(object@package_name) < object@package_version){
					txt <- sprintf('R package %s must have version >= %s', object@package_name, object@package_version)
					return(txt)
				}
			}
		}else if (inherits(object, 'PythonPackage')){

			res <- sprintf('pip show %s', object@package_name) %>% 
				system(intern = TRUE)

			is_available <- is.null(attr(res, 'status'))

			if (!is_available){
				txt <- sprintf('Python package %s is not available', object@package_name)
				return(txt)
			}else{
	
				version <- sprintf('pip show %s', object@package_name) %>% 
					system(intern = TRUE)
				version <- gsub('Version: ', '', version[2])
	
				if (version < object@package_version){
					txt <- sprintf('Python  package %s must have version >= %s', object@package_name, object@package_version)
					return(txt)
				}
			}
		}
	}
	return(NULL)
}

.check_dependences <- function(object){

	res <- lapply(1:length(object@dependences), function(i){
		check_package(object@dependences[[i]])
	})

	res <- res[!sapply(res, is.null)]
	if (length(res) == 0){
		TRUE
	}else{
		for (i in 1:length(res)){
			res[[i]] %>% message()
		}
		return('missing packages')
	}
}


#' RPackage
#' 
#' The base params object for R packages
#'
#' @slot package_name The package name
#' @slot package_version The package version
#'
setClass(
	'RPackage', 
	representation(
    package_name = "character",
		package_version = 'character'
	)
)


#' PythonPackage
#' 
#' The base params object for python packages
#'
#' @slot package_name The package name
#' @slot package_version The package version
#' @slot conda The conda environment
#'
setClass(
	'PythonPackage', 
	representation(
    package_name = "character",
		package_version = 'character'
	)
)

#' BaseMerge 
#' 
#' set base merge params object
#'
#' @slot batch character name of batch in dataset metadata
#' @slot dimreduc_names name of dimension reduction in metadata, i.e. "PCA"
#' @slot return data type to return, can be Seurat or SingleCellExperiment 
#' @slot npcs number of principle components to use in dimension reduction
#' @slot latent whether or not use latent representation to construct the neighboring graph
#'
setClass(
	'BaseMerge', 
	representation(
		batch = "character",
    dimreduc_names = "character",
    return = "character",
		dependences = 'list',
		name = 'character',
		npcs = 'integer',
		latent = 'logical',
		umap_name = 'character',
		umap_key = 'character',
		umap_dim = 'integer',
		snn_name = 'character',
		knn_name = 'character'
	),
  prototype(
		batch = "batch",
		dimreduc_names = c(
			"PCA" = "pca",
			"UMAP" = "umap",
			"tSNE" = "tsne"
		),
		return = "SingleCellExperiment",
		npcs = 20L,
		latent = FALSE,
		umap_dim = 2L
	),
	validity = .check_dependences
)


setMethod('initialize', 'BaseMerge', function(.Object, ...){
	.Object <- callNextMethod()
	.Object@umap_name <- sprintf('%sUMAP', .Object@name)
	.Object@umap_key <- sprintf('%sUMAP_', .Object@name)
	.Object@snn_name <- sprintf('%sSNN', .Object@name)
	.Object@knn_name <- sprintf('%sKNN', .Object@name)
	return(.Object)
})


#' SeuratMerge
#'
#' class for merging dataset with Seurat
#'
#' @slot seed value to set for seed for umap
#' @slot dims number of dimensions used in seurat preprocessing
#' @slot k.weight weight for neighbor function
#'
setClass(
	'SeuratMerge', 
	representation(
    seed = "integer",
    dims = "integer",
    k.weight = "numeric"
	),
	contains = 'BaseMerge',
  prototype(
		seed = 10L,
		dims = 20L,
		k.weight = 100,
		name = 'Seurat'
	)
)



#' BasePreprocess
#' 
#' set base normalize params object
#'
setClass(
	'BasePreprocess', 
	representation(
    min_cells = "integer",
    min_genes = "integer",
		norm_data = "logical",
    scaling = "logical",
    norm_method = "character",
    scale_factor = "numeric"
	),
  prototype(
    min_cells = 10L,
    min_genes = 300L,
		norm_data = TRUE,
		scaling = TRUE,
		norm_method = "LogNormalize",
		scale_factor = 10000
	)
)

#' SeuratPreprocess
#'
#' normalization class for Seurat
#'
#' @slot norm_data boolean check to run data normalization
#' @slot scaling boolean check to scale data
#' @slot regressUMI boolean check to regress the UMI counts
#' @slot min_cells threshold to set the minimum number of cells a gene needs to be present in
#' @slot min_genes threshold of minimum number of genes in cells
#' @slot norm_method set normalization method, default is "LogNormalize"
#' @slot scale_factor factor for scaling
#' @slot vars_to_regress list of names or name of metadata elements to regress
#' @slot numVG number of variable genes for integration
#' @slot numHVG number of highly variable genes to select for downstream analysis
#'
setClass(
	'SeuratPreprocess', 
	representation(
    regressUMI = "logical",
		vars_to_regress = 'character',
    numVG = "integer",
    numHVG = "integer",
		selection.method = 'character'
	),
	contains = c('BasePreprocess'),
  prototype(
		regressUMI = FALSE,
		vars_to_regress = 'nUMI',
		numVG = 300L,
		numHVG = 2000L,
		selection.method = 'vst'
	)
)


#' ScanpyPreprocess
#'
#' @export
setClass(
	'ScanpyPreprocess', 
	representation(
    svd_solver = "character",
    nhvg = "integer",
    n_neighbors = "integer",
		min_mean = 'numeric',
		max_mean = 'integer',
		min_disp = 'numeric'
	),
	contains = 'BasePreprocess',
  prototype(
		min_genes = 200L,
		min_cells  = 3L,
		svd_solver = "arpack",
		nhvg = 2000L,
		n_neighbors = 10L,
		min_mean = 0.0125, 
		max_mean = 3L, 
		min_disp = 0.5
	)
)

#' SeuratParams
#'
#' base class for Seurat processing
#'
#' @slot altExp_names name(s) for additional experiment object in dataset to keep in downstream analysis
#'
#'
#' @export
#'
setClass(
	'SeuratParams', 
	representation(
		altExp_names = "character"
	),
  prototype(
		altExp_names = "RNA",
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0')
		)
	),
	contains = c('SeuratPreprocess', 'SeuratMerge')
)


#' HarmonyMerge
#'
#' merge class for Harmony
#'
#' @slot theta_harmony diversity clustering penalty parameter, larger values increase diversity
#' @slot seed set seed value for umap clustering
#' @slot num_clust number of clusters to use in harmony integration
#' @slot max_iter_cluster maximum number of learning iterations per cluster
#'
#'
#' @export
setClass(
	'HarmonyMerge', 
	representation(
    theta_harmony = "numeric",
    num_clust = "integer",
    max_iter_cluster = "integer",
    seed = "numeric"
	),
	contains = c('BaseMerge'),
  prototype(
		theta_harmony = 2,
		num_clust = 50L,
		max_iter_cluster = 100L,
		seed = 10,
		name = 'Harmony'
	)
)

#' HarmonyParams
#'
#' base parameters class for harmony integration
#'
#'
#' @export
setClass(
	'HarmonyParams', 
	representation(
	),
	prototype(
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0'),
			new('RPackage', package_name = 'harmony', package_version = '0.1.0')
		)
	),
	contains = c('SeuratPreprocess', 'HarmonyMerge')
)

#' UncorrectedMerge
#'
#' class for merging dataset without any correction 
#'
setClass(
	'UncorrectedMerge', 
	representation(
	),
	prototype(
		name = 'Uncorrected'
	),
	contains = c('BaseMerge')
)


#' UncorrectedParams
#'
#' parameters class for uncorrected integration
#'
#' @slot hvg boolean to check if hvg selection should be performed
#'
#'
#' @export
setClass(
	'UncorrectedParams', 
	representation(
	),
	prototype(
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0')
		)
	),
	contains = c('SeuratPreprocess', 'UncorrectedMerge')
)

#' FastMNNMerge
#'
#' class for merging dataset with fastMNN
#'
setClass(
	'FastMNNMerge', 
	representation(
    n_neighbors = "integer"
	),
	contains = c('BaseMerge'),
  prototype(
		n_neighbors = 20L,
		name = 'FastMNN'
	)
)


#' FastMNNParams
#'
#' parameters class for fastMNN integration
#'
#' @export
setClass(
	'FastMNNParams', 
	representation(
	),
	prototype(
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0'),
			new('RPackage', package_name = 'SeuratWrappers', package_version = '0.3.0')
		)
	),
	contains = c('SeuratPreprocess', 'FastMNNMerge')
)


#' LigerMerge
#'
#' parameters class for liger merge
#'
#' @slot nrep number of restarts to perform
#' @slot lambda Regularization parameter. Larger values penalize dataset-specific effects more strongly
#'
#'
#' @export
setClass(
	'LigerMerge', 
	representation(
    nrep = "integer",
    lambda = "numeric"
	),
	contains = c('BaseMerge'),
  prototype(
		nrep = 3L,
		lambda = 5,
		name = 'LIGER'
	)
)

#' LigerParams
#'
#' parameters class for liger operations
#'
#' @export
setClass(
	'LigerParams', 
	representation(
	),
	prototype(
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0'),
			new('RPackage', package_name = 'rliger', package_version = '1.0.0')
		)
	),
	contains = c('SeuratPreprocess', 'LigerMerge')
)



#' BBKNNMerge
#'
#' @export
setClass(
	'BBKNNMerge', 
	representation(
    save_knn = "logical",
    copy = "logical",
    neighbors_within_batch = "numeric",
    approx = "logical",
    trim = "numeric",
    graph_name = "character",
    nn_name = "character",
    ridge_regress = "logical",
    confounder_key = "character"
	),
	contains = 'BaseMerge', 
  prototype(
		save_knn = TRUE,
		copy = TRUE,
		neighbors_within_batch = 5,
		approx = FALSE,
		trim = 50,
		graph_name = "bbknn_graph",
		nn_name = "bbknn",
		ridge_regress = TRUE,
		confounder_key = "leiden",
		name = 'BBKNN'
	)
)

#' BBKNNParams
#'
#' @export
setClass(
	'BBKNNParams', 
	representation(
	),
	prototype(
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0'),
			new('PythonPackage', package_name = 'anndata', package_version = '0.7.8'),
			new('PythonPackage', package_name = 'scanpy', package_version = '1.8'),
			new('PythonPackage', package_name = 'bbknn', package_version = '1.5.1'),
			new('PythonPackage', package_name = 'leidenalg', package_version = '0.8.8'),
			new('RPackage', package_name = 'zellkonverter', package_version = '1.4.0'),
			new('RPackage', package_name = 'basilisk', package_version = '1.6.0')
		)
	),
	contains = c('SeuratPreprocess', 'BBKNNMerge')
)



#' ScanoramaMerge
#'
#' @export
setClass(
	'ScanoramaMerge',
	representation(
	),
	contains = c('BaseMerge'),
  prototype(
		name = "Scanorama"
	)
)

#' ScanoramaParams
#'
#' @export
setClass(
	'ScanoramaParams', 
	representation(
	),
	prototype(
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0'),
			new('PythonPackage', package_name = 'anndata', package_version = '0.7.8'),
			new('PythonPackage', package_name = 'scanpy', package_version = '1.8'),
			new('PythonPackage', package_name = 'scanorama', package_version = '1.7.1'),
			new('RPackage', package_name = 'zellkonverter', package_version = '1.4.0'),
			new('RPackage', package_name = 'basilisk', package_version = '1.6.0')
		)
	),
	contains = c('SeuratPreprocess', 'ScanoramaMerge')
)


#' scVIMerge
#'
#' @export
setClass(
	'scVIMerge',
	representation(
		seed = 'integer'
	),
	contains = c('BaseMerge'),
  prototype(
		seed = 123L,						
		name = "scVI"
	)
)

#' scVIParams 
#'
#' @export
setClass(
	'scVIParams', 
	representation(
	),
	prototype(
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0'),
			new('PythonPackage', package_name = 'anndata', package_version = '0.7.8'),
			new('PythonPackage', package_name = 'scanpy', package_version = '1.8'),
			new('RPackage', package_name = 'zellkonverter', package_version = '1.4.0'),
			new('RPackage', package_name = 'basilisk', package_version = '1.6.0'),
			new('PythonPackage', package_name = 'scvi-tools', package_version = '0.14.5')
		)
	),
	contains = c('SeuratPreprocess', 'scVIMerge')
)

#' ParamsList
#' @export
#'
setClass(
	'ParamsList',
	contains = 'SimpleList'
)

setClass(
	'EnsembleMerge', 
	representation(
		constituent = 'ParamsList'
	),
	prototype(
		name = 'EnsembleMerge'
	),
	contains = 'BaseMerge'
)
	


setClass(
	'EnsembleMergeParams', 
	representation(
	),
	prototype(
		dependences = list(
			new('PythonPackage', package_name = 'umap-learn', package_version = '0.5.2')
		)
	),
	contains = 'EnsembleMerge'
)
	

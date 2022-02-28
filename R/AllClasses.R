#' @import Matrix
#' @importFrom magrittr %>%
#' @import dplyr

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Validity
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#' @importFrom reticulate py_config py_run_string
#' @importFrom utils packageVersion
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
#'
setClass(
	'PythonPackage', 
	representation(
    package_name = "character",
		package_version = 'character'
	)
)



#' The BasePreprocess class
#'
#' @slot min_cells the minimum number of cells that a gene should be expressed.
#' @slot min_genes the minimum number of genes that a cell should express.
#' @slot norm_data whether or not normalizing the data
#' @slot scaling whether or not scaling the data
#' @slot norm_method the normalization method 
#' @slot scale_factor the scaling factor
#' @slot numHVG number of highly variable genes
#' @slot raw_assay the raw assay field in a Seurat object
#' @slot batch character name of batch in dataset metadata
#'
setClass(
	'BasePreprocess', 
	representation(
    min_cells = "integer",
    min_genes = "integer",
		norm_data = "logical",
    scaling = "logical",
    norm_method = "character",
    scale_factor = "numeric",
    numHVG = "integer",
		raw_assay = 'character',
		batch = "character"
	),
  prototype(
    min_cells = 10L,
    min_genes = 0L,
		norm_data = TRUE,
		scaling = TRUE,
		norm_method = "LogNormalize",
		scale_factor = 10000,
		numHVG = 2000L,
		raw_assay = 'RNA',
		batch = 'batch'
	)
)

#' The SeuratPreprocess class
#'
#' @slot selection.method The gene selection method
#' @slot batchwise whether or not performing batchwise data normalization and HVG selection
#'
setClass(
	'SeuratPreprocess', 
	representation(
		selection.method = 'character',
		batchwise = 'logical'
	),
	contains = c('BasePreprocess'),
  prototype(
		selection.method = 'vst',
		batchwise = FALSE
	)
)


#' ScanpyPreprocess
#'
#' @export
#'
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


#' The MethodList class
#'
#' @export
#'
#' @importFrom S4Vectors SimpleList
#'
setClass(
	'MethodList',
	contains = 'SimpleList'
)

#' The SeuratList class
#'
#' @export
#'
#' @importFrom S4Vectors SimpleList
#' @importFrom methods is
#'
setClass(
	'SeuratList',
	contains = 'SimpleList',
	validity = function(object){
		valid <- sapply(object, is, 'Seurat')
		if (any(!valid)){
			sprintf('all elememts must be Seurat objects') %>% message()
			return(FALSE)
		}

		# also need to make sure the the dimensions and other features match

		TRUE
	}
)


setClass(
	'Params', 
	representation(
		preprocess = 'BasePreprocess',
		constituent = 'MethodList'
	),
	prototype(
	)
)


#' BaseMerge 
#' 
#' The BaseMerge class
#'
#' @slot pca_name name of PCA results in metadata
#' @slot dependences a list of dependend packages
#' @slot name The name of the method used to store the dimension reduction results
#' @slot npcs the size of latent dimension (default: 20L)
#' @slot umap_name the name of the UMAP results
#' @slot umap_key the name of the UMAP key
#' @slot umap_dim the name of the UMAP dimension (default: 2L)
#' @slot snn_name the name of SNN results
#' @slot knn_name the name of the KNN results
#' @slot raw_assay the raw assay field in a Seurat object
#' @slot batch character name of batch in dataset metadata
#'
setClass(
	'BaseMerge', 
	representation(
    pca_name = "character",
		dependences = 'list',
		name = 'character',
		reduction_key = 'character',
		npcs = 'integer',
		umap_name = 'character',
		umap_key = 'character',
		umap_dim = 'integer',
		snn_name = 'character',
		knn_name = 'character',
		seed = 'integer',
		raw_assay = 'character',
		batch = "character"
	),
  prototype(
		pca_name = 'pca',
		npcs = 20L,
		umap_dim = 2L,
		seed = 123L,
		raw_assay = 'RNA',
		batch = 'batch'
	)
)

#' @importFrom methods callNextMethod 
#'
setMethod('initialize', 'BaseMerge', function(.Object, ...){
	.check_dependences(.Object)
	.Object@reduction_key <- sprintf('%s_', .Object@name)
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
#' @slot k.weight weight for neighbor function
#'
setClass(
	'SeuratMerge', 
	representation(
    k.weight = "numeric"
	),
	contains = 'BaseMerge',
  prototype(
		k.weight = 100,
		name = 'Seurat',
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0')
		)
	)
)

#' The HarmonyMerge class
#'
#' @slot theta diversity clustering penalty parameter, larger values increase diversity
#' @slot max_iter_cluster maximum number of learning iterations per cluster
#'
#' @export
#'
setClass(
	'HarmonyMerge', 
	representation(
    theta = "numeric",
    max_iter_cluster = "integer"
	),
	contains = c('BaseMerge'),
  prototype(
		theta = 2,
		max_iter_cluster = 20L,
		name = 'Harmony',
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0'),
			new('RPackage', package_name = 'harmony', package_version = '0.1.0')
		)
	)
)

#' The UncorrectedMerge class
#'
#' @export
#'
setClass(
	'UncorrectedMerge', 
	representation(
	),
	prototype(
		name = 'Uncorrected',
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0')
		)
	),
	contains = c('BaseMerge')
)

#' The FastMNNMerge class
#'
#' @slot n_neighbors number of neighbors used in calculating neighboring graph
#'
#' @export
#'
setClass(
	'FastMNNMerge', 
	representation(
    n_neighbors = "integer"
	),
	contains = c('BaseMerge'),
  prototype(
		n_neighbors = 20L,
		name = 'FastMNN',
		dependences = list(
			new('RPackage', package_name = 'batchelor', package_version = '1.10.0'),
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0')
		)
	)
)

#' The LigerMerge class
#'
#' @export
#'
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
		name = 'LIGER',
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0'),
			new('RPackage', package_name = 'rliger', package_version = '1.0.0')
		)
	)
)

#' The BBKNNMerge class
#'
#' @export
#'
setClass(
	'BBKNNMerge', 
	representation(
	),
	contains = 'BaseMerge', 
  prototype(
		name = 'BBKNN',
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0'),
			new('PythonPackage', package_name = 'anndata', package_version = '0.7.8'),
			new('PythonPackage', package_name = 'scanpy', package_version = '1.8'),
			new('PythonPackage', package_name = 'bbknn', package_version = '1.5.1'),
			new('PythonPackage', package_name = 'leidenalg', package_version = '0.8.8'),
			new('RPackage', package_name = 'zellkonverter', package_version = '1.4.0'),
			new('RPackage', package_name = 'basilisk', package_version = '1.6.0'),
			new('PythonPackage', package_name = 'umap-learn', package_version = '0.5.2')
		)
	)
)

#' The ScanoramaMerge class
#'
#' @export
#'
setClass(
	'ScanoramaMerge',
	representation(
	),
	contains = c('BaseMerge'),
  prototype(
		name = "Scanorama",
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0'),
			new('PythonPackage', package_name = 'anndata', package_version = '0.7.8'),
			new('PythonPackage', package_name = 'scanpy', package_version = '1.8'),
			new('PythonPackage', package_name = 'scanorama', package_version = '1.7.1'),
			new('RPackage', package_name = 'zellkonverter', package_version = '1.4.0'),
			new('RPackage', package_name = 'basilisk', package_version = '1.6.0')
		)
	)
)

#' The scVIMerge class
#'
#' @export
#'
setClass(
	'scVIMerge',
	representation(
	),
	contains = c('BaseMerge'),
  prototype(
		name = "scVI",
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0'),
			new('PythonPackage', package_name = 'anndata', package_version = '0.7.8'),
			new('PythonPackage', package_name = 'scanpy', package_version = '1.8'),
			new('RPackage', package_name = 'zellkonverter', package_version = '1.4.0'),
			new('RPackage', package_name = 'basilisk', package_version = '1.6.0'),
			new('PythonPackage', package_name = 'scvi-tools', package_version = '0.14.5')
		)
	)
)

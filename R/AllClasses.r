#' @import Matrix
#' @importFrom magrittr %>%
#' @import dplyr
#' @importFrom methods new is callNextMethod
#' @importFrom SingleCellExperiment colData rowData
#' @importFrom SummarizedExperiment assays
#' @importFrom reticulate py_config py_run_string
#' @importFrom utils packageVersion
#' @importFrom Seurat CreateAssayObject GetAssayData FindClusters FindNeighbors CreateDimReducObject RunUMAP RunPCA FindTransferAnchors IntegrateEmbeddings SplitObject SCTransform SelectIntegrationFeatures NormalizeData FindVariableFeatures ScaleData FindMarkers
#' @importFrom grDevices png dev.off

setClassUnion('characterOrNULL', c("character", "NULL"))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Validity
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

	if (length(object@dependences) > 0){
		res <- lapply(1:length(object@dependences), function(i){
			check_package(object@dependences[[i]])
		})

		res <- res[!sapply(res, is.null)]
		if (length(res) == 0){
			return(TRUE)
		}else{
			for (i in 1:length(res)){
				res[[i]] %>% message()
			}
			return('missing packages')
		}
	}
	return(TRUE)
}

.check_genome <- function(x){
	stopifnot(x %in% c('hg19', 'mm10'))
	return(TRUE)
}


setClass(
	'RPackage', 
	representation(
    package_name = "character",
		package_version = 'character'
	)
)


setClass(
	'PythonPackage', 
	representation(
    package_name = "character",
		package_version = 'character'
	)
)

setClass(
	'BasePreprocess', 
	representation(
    min_cells = "integer",
    min_genes = "integer",
		raw_assay = 'character',
		batchwise = 'logical',
		batch = "character",
		dependences = 'list',
		check_dependencies = 'logical'
	),
	contains = 'VIRTUAL',				 
  prototype(
    min_cells = 10L,
    min_genes = 300L,
		raw_assay = 'RNA',
		check_dependencies = TRUE
	)
)

setMethod('initialize', 'BasePreprocess', function(.Object, check_dependencies = TRUE, ...){
	if (check_dependencies)
		.check_dependences(.Object)
	.Object <- callNextMethod(.Object, check_dependencies = check_dependencies, ...)	
	if (length(.Object@batch) == 0){
		.Object@batchwise <- FALSE
	}else{
		.Object@batchwise <- TRUE
	}
	.Object
})


setClass(
	'BaseNormalize',
	representation(
		assay_name = 'character',	# output assay name
		raw_assay = 'character',	# input assay name
		preprocess = 'BasePreprocess',
		numHVG = "integer",
		dependences = 'list',
		do.scale = 'logical',
		do.center = 'logical',
		check_dependencies = 'logical'
	),
	contains = 'VIRTUAL',				 
	prototype(
		assay_name = 'RNA',
		raw_assay = 'RNA',
		numHVG = 2000L,
		do.scale = TRUE,
		do.center = TRUE,
		check_dependencies = TRUE
	),
	validity = function(object){
		msg <- NULL
		return(msg)
	}
)

setMethod('initialize', 'BaseNormalize', function(.Object, check_dependencies = TRUE, ...){
	if (check_dependencies)
		.check_dependences(.Object)
	callNextMethod(.Object, check_dependencies = check_dependencies, ...)	
})




setClass(
	'BaseDoubletDetect', 
	representation(
		dbr = 'numeric',
		name = 'character',
		normalize = 'BaseNormalize',
		dependences = 'list',
		check_dependencies = 'logical'
	),
	contains = 'VIRTUAL',				 
	prototype(
		dbr = 0.1,
		check_dependencies = TRUE	
	)
)

#'
setMethod('initialize', 'BaseDoubletDetect', function(.Object, check_dependencies = TRUE, ...){
	if (check_dependencies)
		.check_dependences(.Object)
	 callNextMethod(.Object, check_dependencies = check_dependencies, ...)	
})


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
	'BaseEmbed',
	representation(
		name = 'character',
		reduction_key = 'character',
		reduction_name = 'character',
		ndims = 'integer',
		seed = 'integer',
		normalize = 'BaseNormalize',
		dependences = 'list',
		check_dependencies = 'logical'
	),
	contains = 'VIRTUAL',				 
	prototype(
		ndims = 20L,						
		seed = 123L,
		check_dependencies = TRUE
	)
)

setMethod('initialize', 'BaseEmbed', function(.Object, check_dependencies = TRUE, ...){
	if (check_dependencies)
		.check_dependences(.Object)
	.Object@reduction_key <- sprintf('%s_', .Object@name)
	.Object@reduction_name <- .Object@name
	callNextMethod(.Object, check_dependencies = check_dependencies, ...)
})


setClass(
	'BaseMerge', 
	representation(
		nfeatures = 'integer',
		assay_name = 'character'
	),
	contains = c('VIRTUAL', 'BaseEmbed'),
  prototype(
		nfeatures = 2000L
	),
	validity = function(object){
		msg <- NULL
		if (length(object@normalize@preprocess@batch) == 0)
			msg <- 'object@normalize@preprocess@batch must be specified'
		return(msg)
	}
)

setMethod('initialize', 'BaseMerge', function(.Object, check_dependencies = TRUE, ...){
	.Object@assay_name <- sprintf('%sAssay', .Object@name)
	callNextMethod(.Object, check_dependencies = check_dependencies, ...)
})


setClass(
	'BaseCluster',
	representation(
		name = 'character',
		cluster_name = 'character',
		k = 'integer',
		seed = 'integer',
		dependences = 'list',
		check_dependencies = 'logical'
	),
	contains = 'VIRTUAL',
	prototype(
		seed = 123L,
		check_dependencies = TRUE
	)
)

setMethod('initialize', 'BaseCluster', function(.Object, check_dependencies = TRUE, ...){
	if (check_dependencies)
		.check_dependences(.Object)
	.Object@cluster_name <- .Object@name
	callNextMethod(.Object, check_dependencies = check_dependencies, ...)
})



setClass(
	'BaseAnnotate',
	representation(
		normalize= 'BaseNormalize',
		name = 'character',
		annotate_name = 'character',
		dependences = 'list',
		check_dependencies = 'logical',
		genome = 'character'
	),
	contains = 'VIRTUAL',
	prototype(
		check_dependencies = TRUE	,
		genome = 'hg19'
	),
	validity = function(object){
		msg <- NULL
		if (!object@genome %in% c('hg19', 'mm10'))
			msg <- sprintf('unknown genome: %s', object@genome)
		return(msg)
	}
)

setMethod('initialize', 'BaseAnnotate', function(.Object, check_dependencies = TRUE, ...){
	if (check_dependencies)
		.check_dependences(.Object)
	.Object@annotate_name <- .Object@name
	callNextMethod(.Object, check_dependencies = check_dependencies, ...)
})


setClass(
	'BaseGeneMarkers',
	representation(
		check_dependencies = 'logical',
		genome = 'character',
		name = 'character'
	),
	contains = 'VIRTUAL',
	prototype(
		check_dependencies = TRUE,
		genome = 'hg19'
	),
)


setClass(
	'BaseAmbientRNARemoval', 
	representation(
		name = 'character',
		seed = 'integer',
		normalize = 'BaseNormalize',
		dependences = 'list',
		check_dependencies = 'logical'
	),
	contains = 'VIRTUAL',				 
	prototype(
		seed = 123L,
		check_dependencies = TRUE	
	)
)

setMethod('initialize', 'BaseAmbientRNARemoval', function(.Object, check_dependencies = TRUE, ...){
	if (check_dependencies)
		.check_dependences(.Object)
	 callNextMethod(.Object, check_dependencies = check_dependencies, ...)	
})


setClass(
	'BaseReferenceMap',
	representation(
		reduction_key = 'character',
		reduction_name = 'character',
		normalize_query = 'BaseNormalize',
		normalize_atlas = 'BaseNormalize',
		ndims = 'integer',
		dependences = 'list',
		check_dependencies = 'logical',
		name = 'character'
	),
	contains = 'VIRTUAL',
	prototype(
		ndims = 20L,
		check_dependencies = TRUE
	),
)

setMethod('initialize', 'BaseReferenceMap', function(.Object, check_dependencies = TRUE, ...){
	if (check_dependencies)
		.check_dependences(.Object)
	.Object@reduction_key <- sprintf('%s_', .Object@name)
	.Object@reduction_name <- .Object@name
	 callNextMethod(.Object, check_dependencies = check_dependencies, ...)	
})


setClass(
	'BaseDETest',
	representation(
		cluster = 'BaseCluster',
		control = 'characterOrNULL',
		treatment = 'characterOrNULL',
		seed = 'integer',
		dependences = 'list',
		check_dependencies = 'logical',
		name = 'character',
		fc_name = 'character'
	),
	contains = 'VIRTUAL',
	prototype(
		control = NULL,
		treatment = NULL,
		seed = 1L,
		check_dependencies = TRUE
	),
)

setMethod('initialize', 'BaseDETest', function(.Object, check_dependencies = TRUE, ...){
	if (check_dependencies)
		.check_dependences(.Object)
	.Object@fc_name <- .Object@name
	 callNextMethod(.Object, check_dependencies = check_dependencies, ...)	
})


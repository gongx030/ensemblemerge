#' The BaseEmbed class
#'
setClass(
	'BaseEmbed',
	representation(
		name = 'character',
		reduction_key = 'character',
		reduction_name = 'character',
		ndims = 'integer',
		seed = 'integer',
		dependences = 'list',
		check_dependencies = 'logical'
	),
	prototype(
		ndims = 20L,						
		seed = 123L,
		check_dependencies = TRUE
	)
)

#' @importFrom methods callNextMethod 
#'
setMethod('initialize', 'BaseEmbed', function(.Object, check_dependencies = TRUE, ...){
	if (check_dependencies)
		.check_dependences(.Object)
	.Object@reduction_key <- sprintf('%s_', .Object@name)
	.Object@reduction_name <- .Object@name
	callNextMethod(.Object, check_dependencies = check_dependencies, ...)
})


#' The PCAEmbed class
#'
setClass(
	'PCAEmbed',
	representation(
	),
	contains = c('BaseEmbed'),
	prototype(
		name = 'PCAEmbed'
	)
)


#' Embed a Seurat object with PCA
#'
#' @param data a Seurat object
#' @param params a PCAEmbed object
#' @param ... Additional arguments
#' @return returns a data object with PCA embedding
#' @importFrom methods is
#' @importFrom Seurat ScaleData RunPCA
#'
setMethod(
	'Embed',
	signature(
		x = 'Seurat',
		params = 'PCAEmbed'
	),
	function(
		x,
		params,
		...
	){
		# to be implemented
		# 1. params@ndims should not be greater than # cells
		# stopifnot(valid(x, params))	

		x <- ScaleData(object = x, verbose = FALSE)
		x <- RunPCA(
			x,
			npcs = params@ndims,
			reduction.name = params@reduction_name,
			reduction.key = params@reduction_key,
			verbose = FALSE
		)
		x
	}
)

#' Embed a SeuratList
#'
#' @param data a SeuratList object
#' @param params a PCAEmbed object
#' @param ... Additional arguments
#' @return returns a SeuratList object with latent embedding
#'
setMethod(
	'Embed',
	signature(
		x = 'SeuratList',
		params = 'BaseEmbed'
	),
	function(
		x,
		params,
		...
	){

		for (i in 1:length(x)){
			x[[i]] <- Embed(x[[i]], params, ...)
		}
		x
	}
)

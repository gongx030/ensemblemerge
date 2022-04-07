.check_genome <- function(x){

	stopifnot(x %in% c('hg19', 'mm10'))
	return(TRUE)
}

#' The BaseAnnotate class
#'
setClass(
	'BaseAnnotate',
	representation(
		dependences = 'list',
		check_dependencies = 'logical',
		genome = 'character'
	),
	prototype(
		check_dependencies = TRUE	,
		genome = 'hg19'
	)
)

#' @importFrom methods callNextMethod 
#'
setMethod('initialize', 'BaseAnnotate', function(.Object, check_dependencies = TRUE, ...){
	if (check_dependencies)
		.check_dependences(.Object)
	.check_genome(.Object@genome)
	callNextMethod(.Object, check_dependencies = check_dependencies, ...)
})


#' The scCATCHAnnotate class
#'
setClass(
	'scCATCHAnnotate',
	representation(
	),
	contains = c('BaseAnnotate'),
	prototype(
		dependences = list(
			new('RPackage', package_name = 'scCATCH', package_version = '3.0')
		)
	)
)

#' Annotate a Seurat object with scCATCH (https://cran.r-project.org/web/packages/scCATCH/index.html)
#'
#' @param data a Seurat object
#' @param params a scCATCHAnnotate object
#' @param ... Additional arguments
#' @return returns a data object with cell-wise annotation
#' @importFrom methods is
#'
setMethod(
	'Annotate',
	signature(
		x = 'Seurat',
		params = 'scCATCHAnnotate'
	),
	function(
		x,
		params,
		...
	){
		# stopifnot(valid(x, params))
		browser() 	
	}
)

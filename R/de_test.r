#' DE (differential expression) test for a SeuratList object 
#'
#' @param x a SeuratList object
#' @param params a BaseDETest object
#' @param ... Additional arguments
#' @return returns a SeuratList object with DE results
#' @export
#'
setMethod(
	'DETest',
	signature(
		x = 'SeuratList',
		params = 'BaseDETest'
	),
	function(
		x,
		params,
		...
	){

		for (i in 1:length(x)){
			x[[i]] <- DETest(x[[i]], params, ...)
		}
		x
	}
)

setClass(
	'SeuratDETest',
	representation(
		logfc.threshold = 'numeric',
		test.use = 'character',
		min.pct = 'numeric',
		only.pos = 'logical',
		min.cells.feature = 'integer',
		min.cells.group = 'integer'
	),
	contains = 'BaseDETest',
	prototype(
		name = 'SeuratDETest',
		logfc.threshold = 0.25,
		test.use = 'wilcox',
		min.pct = 0.1,
		only.pos = FALSE,
		min.cells.feature = 3L,
		min.cells.group = 3L
	),
	validity = function(object){
		msg <- NULL
		return(msg)	
	}
)

setMethod('initialize', 'SeuratDETest', function(.Object, ...){
	callNextMethod(.Object, ...)
})



#' DE test of a Seurat objects by the default Seurat pipeline (https://satijalab.org/seurat/articles/de_vignette.html)
#'
#' @param x a Seurat object
#' @param params a SeuratDETest object 
#' @param ... Additional arguments
#' @return a data frame with DE testing results
#' @export
#'
setMethod(
	'DETest',
	signature(
		x = 'Seurat',
		params = 'SeuratDETest'
	),
	function(
		x,
		params,
		...
	){

		stopifnot(valid(x, params))

		x@active.assay <- params@cluster@embedding@normalize@assay_name
		sprintf('DETest | set active assay to %s', x@active.assay) %>% message()

		markers <- FindMarkers(
			x, 
			slot = 'data',
			ident.1 = params@control,
			ident.2 = params@treatment,
#			reduction = params@cluster@embedding@reduction_name,
			reduction = NULL,
			features = NULL,
			logfc.threshold = params@logfc.threshold,
			test.use = params@test.use,
			min.pct = params@min.pct,
			min.diff.pct = -Inf,
			verbose = TRUE,
			only.pos = FALSE,
			max.cells.per.ident = Inf,
			random.seed = params@seed,
			latent.vars = NULL,
			min.cells.feature = params@min.cells.feature,
			min.cells.group = params@min.cells.group,
			pseudocount.use = 1,
			mean.fxn = NULL,
			fc_name = params@fc.name,
			base = 2,
			densify = FALSE
		)
		markers
	}
)

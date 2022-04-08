#' The scCATCHAnnotate class
#'
setClass(
	'scCATCHAnnotate',
	representation(
		cluster = 'BaseCluster',
		tissue = 'character',
		species = 'character',
		cancer = 'character',
		cell_min_pct = 'numeric',
		logfc = 'numeric',
		pvalue = 'numeric'

	),
	contains = c('BaseAnnotate'),
	prototype(
		name = 'scCATCHAnnotate',
		cancer = 'Normal',
		dependences = list(
			new('RPackage', package_name = 'scCATCH', package_version = '3.0')
		),
		cell_min_pct = 0.25,
		logfc = 0.25,
		pvalue = 0.05
	),
	validity = function(object){
		msg <- NULL
		if (!object@genome %in% c('hg19', 'mm10'))
			msg <- sprintf('unknown genome: %s', object@genome)
		return(msg)
	}
)

#' @importFrom methods callNextMethod 
#'
setMethod('initialize', 'scCATCHAnnotate', function(.Object, check_dependencies = TRUE, ...){

	.Object <- callNextMethod(.Object, check_dependencies = check_dependencies, ...)

	if (.Object@genome == 'mm10'){
		.Object@species <- 'Mouse'
	}else if (.Object@genome == 'hg19'){
		.Object@species <- 'Human'
	}
	.Object
})


#' Annotate a Seurat object with scCATCH (https://cran.r-project.org/web/packages/scCATCH/index.html)
#'
#' @param x a Seurat object
#' @param params a scCATCHAnnotate object
#' @param ... Additional arguments
#' @return returns a data object with cell-wise annotation
#' @importFrom rlang .data
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

		# yet to be implemented
		# stopifnot(valid(x, params))

		tissues <- scCATCH::cellmatch %>%
			filter(.data$cancer == params@cancer & .data$species == params@species) %>%
			pull(.data$tissue) %>%
			unique()

		if (length(params@tissue) == 0)
			params@tissue <- tissues

		if (!all(params@tissue %in% tissues))
			stop(sprintf('params@tissue must be one of the following: %s', paste(tissues, collapse = ',')))

		raw_assay <- params@cluster@embedding@preprocess@raw_assay
		cls <- x[[params@cluster@cluster_name]][, 1] %>%
			as.character()

		obj <- scCATCH::createscCATCH(data = x@assays[[raw_assay]]@data, cluster = cls)
		obj <- scCATCH::findmarkergene(
			object = obj, 
			species = params@species, 
			marker = scCATCH::cellmatch, 
			tissue = params@tissue, 
			cell_min_pct = params@cell_min_pct,
			logfc = params@logfc,
			pvalue = params@pvalue,
			cancer = params@cancer,
			verbose = FALSE
		)
		obj <- scCATCH::findcelltype(object = obj, verbose = FALSE)

		rownames(obj@celltype) <- obj@celltype[, 'cluster']
		x@meta.data[[params@annotate_name]] <- obj@celltype[obj@meta$cluster, 'cell_type']
		x
	}
)


#' Annotate a SeuratList
#'
#' @param x a SeuratList object
#' @param params a BaseAnnotate object
#' @param ... Additional arguments
#' @return returns a SeuratList object with cell-wise annotation
#'
setMethod(
	'Annotate',
	signature(
		x = 'SeuratList',
		params = 'BaseAnnotate'
	),
	function(
		x,
		params,
		...
	){

		for (i in 1:length(x)){
			x[[i]] <- Annotate(x[[i]], params, ...)
		}
		x
	}
)

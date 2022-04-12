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

#' The SCINAAnnotate class
#'
setClass(
	'SCINAAnnotate',
	representation(
		species = 'character',
		cancer = 'character',
		max_iter = 'integer',
		convergence_n = 'integer',
		convergence_rate = 'numeric',
		sensitivity_cutoff = 'numeric',
		rm_overlap = 'logical',
		allow_unknown = 'logical'
	),
	contains = c('BaseAnnotate'),
	prototype(
		name = 'SCINAAnnotate',
		cancer = 'Normal',
		max_iter = 100L,
		convergence_n = 10L,
		convergence_rate = 0.999,
		sensitivity_cutoff = 0.9,
		rm_overlap = FALSE,
		allow_unknown = TRUE,
		dependences = list(
			new('RPackage', package_name = 'SCINA', package_version = '1.2.0'),
			new('RPackage', package_name = 'scCATCH', package_version = '3.0')
		)
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
setMethod('initialize', 'SCINAAnnotate', function(.Object, check_dependencies = TRUE, ...){
	.Object <- callNextMethod(.Object, check_dependencies = check_dependencies, ...)
	if (.Object@genome == 'mm10'){
		.Object@species <- 'Mouse'
	}else if (.Object@genome == 'hg19'){
		.Object@species <- 'Human'
	}
	.Object
})


#' Annotate a Seurat object with SCINA (https://github.com/jcao89757/SCINA)
#'
#' @param x a Seurat object
#' @param params a SCINAAnnotate object
#' @param ... Additional arguments
#' @return returns a data object with cell-wise annotation
#' @importFrom rlang .data
#'
setMethod(
	'Annotate',
	signature(
		x = 'Seurat',
		params = 'SCINAAnnotate'
	),
	function(
		x,
		params,
		...
	){

		# yet to be implemented
		# stopifnot(valid(x, params))

		signatures <- scCATCH::cellmatch %>%
			filter(.data$cancer == params@cancer & .data$species == params@species) %>%
			select(.data$gene, .data$celltype)

		signatures <- split(signatures$gene, list(signatures$celltype))

		raw_assay <- params@preprocess@raw_assay

		results <- SCINA::SCINA(
			x@assays[[raw_assay]]@data, 
			signatures,
			max_iter = params@max_iter,
			convergence_n = params@convergence_n,
			convergence_rate = params@convergence_rate,
			sensitivity_cutoff = params@sensitivity_cutoff,
			rm_overlap = params@rm_overlap,
			allow_unknown = params@allow_unknown
		)

		x@meta.data[[params@annotate_name]] <- results$cell_labels
		x
	}
)

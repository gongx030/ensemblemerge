#' Preprocess a SeuratList objects by the Seurat pipeline
#'
#' @param x a SeuratList object
#' @param params a BasePreprocess object 
#' @param ... Additional arguments
#' @return a SeuratList object 
#' @export
setMethod(
	'Preprocess',
	signature(
		x = 'SeuratList',
		params = 'BasePreprocess'
	),
	function(
		x,
		params,
		...
	){
		for (i in 1:length(x)){
			x[[i]] <- Preprocess(x[[i]], params, ...)
		}
		x	
	}
)

setClass(
	'SeuratPreprocess', 
	representation(
	),
	contains = c('BasePreprocess'),
  prototype(
	)
)

#' @importFrom methods callNextMethod
#'
setMethod('initialize', 'SeuratPreprocess', function(.Object, check_dependencies = TRUE, ...){
	if (check_dependencies)
		.check_dependences(.Object)
	.Object <- callNextMethod(.Object, check_dependencies = check_dependencies, ...)
	.Object
})


#' Preprocess a Seurat objects 
#'
#' @param x a Seurat object
#' @param params a BasePreprocess object 
#' @param ... Additional arguments
#' @return a Seurat object (if there is only one batch), or a SeuratList
#' @export
#' @importFrom Seurat GetAssayData
#'
setMethod(
	'Preprocess',
	signature(
		x = 'Seurat',
		params = 'BasePreprocess'
	),
	function(
		x,
		params,
		...
	){

		stopifnot(valid(x, params))

		rc <- rowSums(GetAssayData(x, 'data') > 0)
		cc <- colSums(GetAssayData(x, 'data') > 0)

		invalid <- rc < params@min_cells
		if (any(invalid)){
			x <- x[!invalid, ]
			sprintf('Preprocess | removing %d genes that are expressed in <%d (min_cells) cells', sum(invalid), params@min_cells) %>% message()
		}

		invalid <- cc < params@min_genes
		if (any(invalid)){
			x <- x[, !invalid]
			sprintf('Preprocess | removing %d cells that have <%d (min_genes) detected genes', sum(invalid), params@min_genes) %>% message()
		}

		x	
	}
)


#' Preprocess a SummarizedExperiment objects by the Seurat pipeline
#'
#' Adopted from https://satijalab.org/seurat/articles/integration_introduction.html
#'
#' @param x a Seurat object
#' @param params a SeuratPreprocess object 
#' @param counts the assay field for raw counts in a SingleCellExperiment object (default: 'counts')
#' @param ... Additional arguments
#' @return a Seurat object (if there is only one batch), or a SeuratList
#' @export
#' @importFrom Seurat SplitObject NormalizeData FindVariableFeatures ScaleData SelectIntegrationFeatures CreateSeuratObject
#' @importFrom SingleCellExperiment colData
#'
setMethod(
	'Preprocess',
	signature(
		x = 'SummarizedExperiment',
		params = 'SeuratPreprocess'
	),
	function(
		x,
		params,
		counts = 'counts',
		...
	){

		# This need to be implemented
		# stopifnot(valid(data, params))

		seurat <- CreateSeuratObject(
			counts = assays(x)[[counts]], # Unnormalized data such as raw counts or TPMs
			meta.data = as.data.frame(colData(x)),
			assay = params@raw_assay
		)

		if (ncol(rowData(x)) > 0){
			d <- as.data.frame(rowData(x))
			seurat[[params@raw_assay]]@meta.features <- d
		}


		Preprocess(seurat, params, ...)

	}
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



#' Preprocess a Seurat object by the Scanpy pipeline
#'
#' @param x a Seurat object
#' @param params a ScanpyPreprocess object
#' @param merge a BaseMerge object
#' @param ... Additional arguments
#' @return returns a Seurat object
#'
#setMethod(
#	'Preprocess',
#	signature(
#		x = 'Seurat',
#		params = 'ScanpyPreprocess',
#		merge = 'BaseMerge'
#	),
#	function(
#		x,
#		params,
#		merge,
#		...
#	){
#
#		x <- as.SingleCellExperiment(x)
#		Preprocess(params, x, ...)
#	}
#)

#' Preprocess a SingleCellExperiment object by the Scanpy pipeline
#'
#' @param x a SingleCellExperiment object
#' @param params a ScanpyPreprocess object
#' @return returns a SingleCellExperiment object
#' @param ... Additional arguments
#' @importFrom reticulate import
#'
#setMethod(
#	'Preprocess',
#	signature(
#		x = 'SingleCellExperiment',
#		params = 'ScanpyPreprocess'
#	),
#	function(
#		x,
#		params,
#		...
#	){
#
#		h5ad_file <- tempfile(fileext = '.h5ad')
#		x = sceasy::convertFormat(x, from = "sce", to = "anndata", out = h5ad_file)
#		sc <- import("scanpy")
#		sc$pp$filter_cells(x, min_genes = params@min_genes)
#		sc$pp$filter_genes(x, min_cells = params@min_cells)
#		sc$pp$normalize_total(x, target_sum = params@scale_factor)
#		sc$pp$log1p(x)
#		sc$pp$highly_variable_genes(x, min_mean = params@min_mean, max_mean= params@max_mean, min_disp = params@min_disp)
#		sc$pp$scale(x)
#		sc$tl$pca(x)
#		x$write(filename = h5ad_file)
#		x <- zellkonverter::readH5AD(h5ad_file, reader = 'R')
#		rowData(x)$highly_variable <- as.logical(rowData(x)$highly_variable)
#		unlink(h5ad_file)
#		x	
#	}
#)

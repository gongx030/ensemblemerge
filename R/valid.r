#' Validate a Seurat object before preprocessing
#'
#' @param x a Seurat object
#' @param params a BasePreprocess object
#' @param ... Additional arguments

#' @importFrom Seurat DefaultAssay<-
#' @importFrom Matrix rowSums colSums
#
setMethod(
	'valid',
	signature(
		x = 'Seurat',
		params = 'BasePreprocess'
	),
	function(
		x,
		params,
		...
	){

		if (any(duplicated(colnames(x))))
			stop('there are duplicated column names in the input data')

		if (!params@raw_assay %in% names(x@assays))
			stop(sprintf('The raw counts are not present in x@assays$%s', params@raw_assay))

		DefaultAssay(x) <- params@raw_assay

		if (!params@batch %in% colnames(x@meta.data))
			stop(sprintf('the batch field is not present in x[["%s"]]', params@batch))

		return(TRUE)

	}
)

#' Validate a Seurat object before merging
#'
#' @param x a Seurat object
#' @param params a BaseMerge object
#' @param ... Additional arguments
#'
#' @importFrom Seurat DefaultAssay<-
#' @importFrom Matrix rowSums colSums
#
setMethod(
	'valid',
	signature(
		x = 'Seurat',
		params = 'BaseMerge'
	),
	function(
		x,
		params,
		...
	){

		if (any(duplicated(colnames(x))))
			stop('there are duplicated column names in the input data')

		if (!params@raw_assay %in% names(x@assays))
			stop(sprintf('The raw counts are not present in x@assays$%s', params@raw_assay))

		DefaultAssay(x) <- params@raw_assay

		if (!params@batch %in% colnames(x@meta.data))
			stop(sprintf('the batch field is not present in x[["%s"]]', params@batch))

		return(TRUE)

	}
)

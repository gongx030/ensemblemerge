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

		if (any(duplicated(rownames(x))))
			stop('there are duplicated row names in the input data')

		if (!params@raw_assay %in% names(x@assays))
			stop(sprintf('The raw counts are not present in x@assays$%s', params@raw_assay))

		DefaultAssay(x) <- params@raw_assay

		if (length(params@batch) > 0){
			if (!params@batch %in% colnames(x@meta.data))
				stop(sprintf('the batch field is not present in x[["%s"]]', params@batch))
		}

		return(TRUE)

	}
)

#' Validate a Seurat object before normalization
#'
#' @param x a Seurat object
#' @param params a BaseNormalize object
#' @param ... Additional arguments
#
setMethod(
	'valid',
	signature(
		x = 'Seurat',
		params = 'BaseNormalize'
	),
	function(
		x,
		params,
		...
	){

		valid(x, params@preprocess)
		return(TRUE)
	}
)

#' Validate a SeuratList object before merging
#'
#' @param x a SeuratList object
#' @param params a BaseMerge object
#' @param ... Additional arguments
#
setMethod(
	'valid',
	signature(
		x = 'SeuratList',
		params = 'BaseMerge'
	),
	function(
		x,
		params,
		...
	){

		for (i in 1:length(x)){
			valid(x[[i]], params)
		}
		return(TRUE)
	}
)

#' Validate a Seurat object before merging
#'
#' @param x a Seurat object
#' @param params a BaseMerge object
#' @param ... Additional arguments
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

		valid(x, params@normalize)
		return(TRUE)
	}
)

#' Validate a Seurat object before DE test
#'
#' @param x a Seurat object
#' @param params a BaseDETest object
#' @param ... Additional arguments
#
setMethod(
	'valid',
	signature(
		x = 'Seurat',
		params = 'BaseDETest'
	),
	function(
		x,
		params,
		...
	){
		if (!is.null(params@treatment)){
			is_valid <- params@treatment %in% x[[params@cluster@cluster_name]][, 1]
			if (any(!is_valid)){
				sprintf('the following in params@treatment is not present in x[[%s]]: %s', params@cluster@cluster_name, paste(params@treatment[!is_valid], collapse = ',')) %>% message()
			}
		}

		if (!is.null(params@control)){
			is_valid <- params@control %in% x[[params@cluster@cluster_name]][, 1]
			if (any(!is_valid)){
				sprintf('the following in params@control is not present in x[[%s]]: %s', params@cluster@cluster_name, paste(params@control[!is_valid], collapse = ',')) %>% message()
			}
		}

		return(TRUE)
	}
)


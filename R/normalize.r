#' Normalize a SeuratList object 
#'
#' @param x a SeuratList object
#' @param params a BaseNormalize object
#' @param ... Additional arguments
#' @return returns a SeuratList object with normalized read counts
#' @export
#'
setMethod(
	'Normalize',
	signature(
		x = 'SeuratList',
		params = 'BaseNormalize'
	),
	function(
		x,
		params,
		...
	){

		for (i in 1:length(x)){
			x[[i]] <- Normalize(x[[i]], params, ...)
		}
		x
	}
)


setClass(
	'SeuratNormalize',
	representation(
		selection.method = 'character',
		norm_method = 'character',
    scaling = "logical",
    scale_factor = "numeric",
		feature_field = 'character',
		output = 'character'
	),
	contains = 'BaseNormalize',
	prototype(
		selection.method = 'vst',
		norm_method = "LogNormalize",
		scaling = TRUE,
		scale_factor = 10000,
		output = 'Seurat'
	),
	validity = function(object){
		msg <- NULL
		if (!object@selection.method %in% c('vst', 'mean.var.plot', 'dispersion'))
			msg <- sprintf('unknown selection.method: %s', object@selection.method)
		return(msg)	
	}
)

setMethod('initialize', 'SeuratNormalize', function(.Object, ...){

	if (.Object@selection.method == 'vst')
		.Object@feature_field <- 'vst.variable'
	else if (.Object@selection.method == 'mean.var.plot')
		.Object@feature_field <- 'mvp.variable'
	else if (.Object@selection.method == 'dispersion')
		.Object@feature_field <- 'mvp.variable'

	.Object@raw_assay <- .Object@preprocess@raw_assay
	.Object@assay_name <- .Object@raw_assay

	callNextMethod(.Object, ...)

})


#' Normalize a Seurat objects by the default Seurat pipeline
#'
#' Adopted from https://satijalab.org/seurat/articles/integration_introduction.html
#'
#' @param x a Seurat object
#' @param params a SeuratNormalize object 
#' @param ... Additional arguments
#' @return a Seurat object (if there is only one batch), or a SeuratList
#' @export
#' @importFrom Seurat SplitObject NormalizeData FindVariableFeatures ScaleData SelectIntegrationFeatures
#'
setMethod(
	'Normalize',
	signature(
		x = 'Seurat',
		params = 'SeuratNormalize'
	),
	function(
		x,
		params,
		...
	){

		x@active.assay <- params@raw_assay

		if (params@preprocess@batchwise){
			x <- SplitObject(x, split.by = params@preprocess@batch)
		}else{
			x <- list(x)
		}

		for (i in 1:length(x)){

			sprintf('Normalize | input assay=%s | output assay=%s', params@raw_assay, params@assay_name) %>% message()

			x[[i]] <- NormalizeData(
				x[[i]],
				normalization.method = params@norm_method,
				scale.factor = params@scale_factor,
				verbose = FALSE
			)

	 	  x[[i]] <- FindVariableFeatures(
				x[[i]],
				selection.method = params@selection.method,
				nfeatures = params@numHVG,
				verbose = FALSE
			)

			x[[i]] <- ScaleData(
				x[[i]],
				vars.to.regress = NULL,
				do.scale = params@do.scale,
				do.center = params@do.center,
				use.umi = FALSE,
				verbose = FALSE
			)
		}

		if (length(x) == 1){
			x[[1]]
		}else{
			new('SeuratList', x)
		}

	}
)

setClass(
	'SCTransformNormalize',
	representation(
		do.correct.umi = 'logical',
		ncells = 'integer'
	),
	contains = 'BaseNormalize',
	prototype(
		assay_name = 'SCT',
		do.correct.umi = TRUE,
		ncells = 2000L
	),
	validity = function(object){
		msg <- NULL
		return(msg)	
	}
)

#' @importFrom methods callNextMethod
#'
setMethod('initialize', 'SCTransformNormalize', function(.Object, ...){
	.Object@raw_assay <- .Object@preprocess@raw_assay
	callNextMethod(.Object, ...)
})


#' Normalize a Seurat objects by SCTransform normalization
#'
#' Adopted from https://satijalab.org/seurat/articles/sctransform_vignette.html
#'
#' @param x a Seurat object
#' @param params a SCTransformNormalize object 
#' @param ... Additional arguments
#' @return a Seurat object (if there is only one batch), or a SeuratList
#' @export
#' @importFrom Seurat SplitObject SCTransform SelectIntegrationFeatures
#' @references Hafemeister, C., Satija, R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome Biol 20, 296 (2019). https://doi.org/10.1186/s13059-019-1874-1
#'
setMethod(
	'Normalize',
	signature(
		x = 'Seurat',
		params = 'SCTransformNormalize'
	),
	function(
		x,
		params,
		...
	){

		if (params@preprocess@batchwise){
			x <- SplitObject(x, split.by = params@preprocess@batch)
		}else{
			x <- list(x)
		}

		for (i in 1:length(x)){

			sprintf('Normalize | input assay=%s | output assay=%s', params@raw_assay, params@assay_name) %>% message()

			x[[i]] <- SCTransform(
				x[[i]],
				assay = params@raw_assay,
				new.assay.name = params@assay_name,
				reference.SCT.model = NULL,
				do.correct.umi = params@do.correct.umi,
				ncells = params@ncells,
				residual.features = NULL,
				variable.features.n = params@numHVG,
				variable.features.rv.th = NULL,
				vars.to.regress = NULL,
				do.scale = params@do.scale,
				do.center = params@do.center,
				conserve.memory = FALSE,
				return.only.var.genes = TRUE,
				verbose = FALSE 
			)
			sprintf('Normalize | set active assay to "%s"', params@assay_name) %>% message()
		}

		if (length(x) == 1){
			x[[1]]
		}else{
			new('SeuratList', x)
		}
	}
)


setClass(
	'scranNormalize',
	representation(
		min_size = 'integer',
		clustering_method = 'character',
		ndims = 'integer',
		graph_fun = 'character'
	),
	contains = 'BaseNormalize',
	prototype(
		min_size = 100L,
		clustering_method = 'igraph',
		ndims = 50L,
		graph_fun = 'walktrap',
		dependences = list(
			new('RPackage', package_name = 'scran', package_version = '1.22.1'),
			new('RPackage', package_name = 'scuttle', package_version = '1.4.0')
		)
	),
	validity = function(object){
		msg <- NULL
		return(msg)	
	}
)

setMethod('initialize', 'scranNormalize', function(.Object, ...){
	.Object@raw_assay <- .Object@preprocess@raw_assay
	callNextMethod(.Object, ...)
})


#' Normalize a Seurat objects by scran's devolution normalization (https://bioconductor.org/packages/devel/bioc/vignettes/scran/inst/doc/scran.html)
#'
#' @param x a Seurat object
#' @param params a scranNormalize object 
#' @param ... Additional arguments
#' @return a Seurat object (if there is only one batch), or a SeuratList
#' @export
#' @references L. Lun, A.T., Bach, K. & Marioni, J.C. Pooling across cells to normalize single-cell RNA sequencing data with many zero counts. Genome Biol 17, 75 (2016). https://doi.org/10.1186/s13059-016-0947-7
#'
setMethod(
	'Normalize',
	signature(
		x = 'Seurat',
		params = 'scranNormalize'
	),
	function(
		x,
		params,
		...
	){


		if (params@preprocess@batchwise){
			x <- SplitObject(x, split.by = params@preprocess@batch)
		}else{
			x <- list(x)
		}

		for (i in 1:length(x)){

			sprintf('Normalize | input assay=%s | output assay=%s', params@raw_assay, params@assay_name) %>% message()

			x[[i]]@active.assay <- params@assay_name

			sce <- SingleCellExperiment(
				assays = list(
					counts = GetAssayData(x[[i]], 'counts')
				),
			)

			clusters <- scran::quickCluster(
				sce,
				min.size = params@min_size,
				method = params@clustering_method,
				use.ranks = FALSE,
				d = params@ndims,
				subset.row = NULL,
				min.mean = NULL,
				graph.fun = params@graph_fun,
				block = NULL
			)

			sce <- scran::computeSumFactors(sce, clusters = clusters)
			sce <- scuttle::logNormCounts(sce)
			dec <- scran::modelGeneVar(sce)
			features <- scran::getTopHVGs(
				dec, 
				n = params@numHVG,
				prop = NULL,
				var.threshold = 0,
				fdr.field = "FDR",
				fdr.threshold = NULL
			)

			x[[i]]@assays[[x[[i]]@active.assay]]@var.features <- features

			x[[i]] <- ScaleData(
				x[[i]],
				vars.to.regress = NULL,
				do.scale = params@do.scale,
				do.center = params@do.center,
				use.umi = FALSE,
				verbose = FALSE
			)
		}

		if (length(x) == 1){
			x[[1]]
		}else{
			new('SeuratList', x)
		}
	}
)

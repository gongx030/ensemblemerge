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

#' @importFrom methods callNextMethod
#'
setMethod('initialize', 'SeuratNormalize', function(.Object, check_dependencies = TRUE, ...){

	if (.Object@selection.method == 'vst')
		.Object@feature_field <- 'vst.variable'
	else if (.Object@selection.method == 'mean.var.plot')
		.Object@feature_field <- 'mvp.variable'
	else if (.Object@selection.method == 'dispersion')
		.Object@feature_field <- 'mvp.variable'

	.Object@assay_name <- .Object@preprocess@raw_assay

	callNextMethod(.Object, check_dependencies = check_dependencies, ...)

})


#' Normalize a Seurat objects by the default Seurat pipeline
#'
#' Adopted from https://satijalab.org/seurat/articles/integration_introduction.html
#'
#' @param x a Seurat object
#' @param params a SeuratNormalize object 
#' @param ... Additional arguments
#' @return a Seurat object (if there is only one batch), or a SeuratList

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

		x@active.assay <- params@assay_name

		if (params@batchwise){
			batch_list <- SplitObject(x, split.by = params@preprocess@batch)
		}else{
			batch_list <- list(x)
		}

		for (i in 1:length(batch_list)){
			batch_list[[i]] <- NormalizeData(
				batch_list[[i]],	
				normalization.method = params@norm_method,
				scale.factor = params@scale_factor,
				verbose = FALSE
			)
 	   	batch_list[[i]]  <- FindVariableFeatures(
				batch_list[[i]],
				selection.method = params@selection.method,
				nfeatures = params@numHVG,
				verbose = FALSE
			)
		}

		if (params@output == 'SeuratList'){

			for (i in 1:length(batch_list)){
				batch_list[[i]] <- ScaleData(
					object = batch_list[[i]], 
					do.scale = params@do.scale,
					do.center = params@do.center,
					use.umi = FALSE,
					verbose = FALSE
				)
			}
			new('SeuratList', batch_list)

		}else if (params@output == 'Seurat'){
			# select features that are repeatedly variable across datasets for integration
			features <- SelectIntegrationFeatures(
				object.list = batch_list,
				nfeatures = params@numHVG
			)
			x@assays[[x@active.assay]]@var.features <- features
			x <- ScaleData(
				object = x, 
				do.scale = params@do.scale,
				do.center = params@do.center,
				use.umi = FALSE,
				verbose = FALSE
			)
			x	
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
setMethod('initialize', 'SCTransformNormalize', function(.Object, check_dependencies = TRUE, ...){
	callNextMethod(.Object, check_dependencies = check_dependencies, ...)
})


#' Normalize a Seurat objects by SCTransform normalization
#'
#' Adopted from https://satijalab.org/seurat/articles/sctransform_vignette.html
#'
#' @param x a Seurat object
#' @param params a SCTransformNormalize object 
#' @param ... Additional arguments
#' @return a Seurat object (if there is only one batch), or a SeuratList

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

		if (params@batchwise){
			vars.to.regress <- params@preprocess@batch
		}else{
			vars.to.regress <- NULL
		}

		x <- SCTransform(
			x,
			assay = params@preprocess@raw_assay,
			new.assay.name = params@assay_name,
			reference.SCT.model = NULL,
			do.correct.umi = params@do.correct.umi,
			ncells = params@ncells,
			residual.features = NULL,
			variable.features.n = params@numHVG,
			variable.features.rv.th = NULL,
			vars.to.regress = vars.to.regress,
			do.scale = params@do.scale,
			do.center = params@do.center,
			conserve.memory = FALSE,
			return.only.var.genes = TRUE,
			verbose = FALSE 
		)

		if (params@output == 'SeuratList'){
			batch_list <- SplitObject(x, split.by = params@preprocess@batch)
			new('SeuratList', batch_list)
		}else if (params@output == 'Seurat'){
			x
		}
	}
)

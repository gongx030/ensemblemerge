
setClass(
	'SeuratMerge', 
	representation(
    k.weight = "numeric",
		normalization.method = 'character',
		reduction = 'character'
	),
	contains = 'BaseMerge',
  prototype(
		name = 'SeuratMerge',
		k.weight = 100,
		reduction = 'cca',
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0')
		)
	)
)

setMethod('initialize', 'SeuratMerge', function(.Object, ...){
	.Object <- callNextMethod(.Object, ...)
	if (is(.Object@normalize, 'SeuratNormalize'))
		.Object@normalization.method <- 'LogNormalize'
	else if (is(.Object@normalize, 'SCTransformNormalize'))
		.Object@normalization.method <- 'SCT'
	.Object
})


#' Batch correction by Seurat integration pipeline (https://satijalab.org/seurat/articles/integration_introduction.html)
#'
#' @param x a SeuratList object 
#' @param params a SeuratMerge object 
#' @param ... Additional arguments
#' @return returns a data object with integrated data
#' @export
#' @importFrom methods is
#' @importFrom Seurat SplitObject SelectIntegrationFeatures FindIntegrationAnchors PrepSCTIntegration IntegrateData
#'
setMethod(
	'Merge',
	signature(
		x = 'SeuratList',
		params = 'SeuratMerge'
	),
	function(
		x,
		params,
		...
	){

		stopifnot(valid(x, params))

		for (i in 1:length(x)){
			x[[i]]@active.assay <- params@normalize@assay_name
		}

		features <- SelectIntegrationFeatures(
			object.list = x,
			nfeatures = params@nfeatures,
			verbose = FALSE
		)

		# when k.weight is greater than number of cell in any batches, it will result such error:
		# IntegrateData error: Error in idx[i, ] <- res[[i]][[1]]
		# https://github.com/satijalab/seurat/issues/3930
		n <- sapply(x, ncol) # number of cells per cluster
		k_weight <- min(n, params@k.weight)

		if (params@normalization.method == 'SCT'){

			x <- PrepSCTIntegration(
				x,
				anchor.features = features,
				verbose = FALSE
			)
		}

		x <- FindIntegrationAnchors(
			object.list = x,
			anchor.features = features,
			reference = NULL,
			scale = TRUE,
			normalization.method = params@normalization.method,
			sct.clip.range = NULL,
			reduction = params@reduction,
			l2.norm = TRUE,
			dims = 1:params@ndims,
			k.anchor = 5,
			k.filter = 200,
			k.score = 30,
			max.features = 200,
			nn.method = "annoy",
			n.trees = 50,
			eps = 0,
			verbose = FALSE 
		)

	 	x <- IntegrateData(
			anchorset = x, 
			new.assay.name = params@assay_name,
			dims = 1:params@ndims,
			k.weight = k_weight,
			normalization.method = params@normalization.method,
			features = NULL,
			features.to.integrate = NULL,
			weight.reduction = NULL,
			sd.weight = 1,
			sample.tree = NULL,
			preserve.order = FALSE,
			eps = 0,
			verbose = TRUE
		)

		x <- ScaleData(
			object = x, 
			vars.to.regress = NULL,
			do.scale = params@normalize@do.scale,
			do.center = params@normalize@do.center,
			use.umi = FALSE,
			verbose = FALSE
		)

	  x <- RunPCA(
			x, 
			npcs = params@ndims, 
			reduction.name = params@reduction_name,,
			reduction.key = params@reduction_key,
			verbose = FALSE
		)

		x@assays[[params@normalize@assay_name]]@var.features <- features

		x
	}
)


setClass(
	'HarmonyMerge', 
	representation(
    theta = "numeric",
    max_iter_cluster = "integer"
	),
	contains = c('BaseMerge'),
  prototype(
		theta = 2,
		max_iter_cluster = 20L,
		name = 'HarmonyMerge',
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0'),
			new('RPackage', package_name = 'harmony', package_version = '0.1.0')
		)
	)
)

#' Merge Seurat objects
#'
#' @param x Seurat object containing single cell counts matrix
#' @param params a HarmonyMerge object 
#' @param ... Additional arguments
#' @return returns a Seurat object with integrated data
#' @export
#'
setMethod(
	'Merge',
	signature(
		x = 'SeuratList',
		params = 'HarmonyMerge'
	),
	function(
		x,
		params,
		...
	){


		stopifnot(valid(x, params))

		for (i in 1:length(x)){
			x[[i]]@active.assay <- params@normalize@assay_name
		}

		features <- SelectIntegrationFeatures(
			object.list = x,
			nfeatures = params@nfeatures,
			verbose = FALSE
		)

		x <- Reduce('merge', x)

		x@assays[[params@normalize@assay_name]]@var.features <- features

		x <- ScaleData(
			object = x, 
			vars.to.regress = NULL,
			do.scale = params@normalize@do.scale,
			do.center = params@normalize@do.center,
			use.umi = FALSE,
			verbose = FALSE
		)

		params_pca <- new('PCAEmbed', normalize = params@normalize, ndims = params@ndims)

		x <- Embed(x, params_pca)

		x <- harmony::RunHarmony(
			object = x, 
			group.by.vars = params@normalize@preprocess@batch, 
			reduction = params_pca@reduction_name,
			theta = params@theta, 
			plot_convergence = FALSE, 
			max.iter.cluster = params@max_iter_cluster,
			assay.use = params@normalize@assay_name,
			reduction.save = params@reduction_name,
			verbose = FALSE 
		)

		x	
	}
)


setClass(
	'FastMNNMerge', 
	representation(
    n_neighbors = "integer"
	),
	contains = c('BaseMerge'),
  prototype(
		n_neighbors = 20L,
		name = 'FastMNN',
		dependences = list(
			new('RPackage', package_name = 'batchelor', package_version = '1.10.0'),
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0')
		)
	)
)

#' Merge Seurat objects
#'
#' @param x Seurat object containing single cell counts matrix
#' @param params a FastMNNMerge object 
#' @param ... Additional arguments
#' @return returns a Seurat object with integrated data
#' @importFrom Seurat VariableFeatures SplitObject RunPCA as.SingleCellExperiment CreateDimReducObject
#' @importFrom SingleCellExperiment reducedDim
#' @export
#'
setMethod(
	'Merge',
	signature(
		x = 'SeuratList',
		params = 'FastMNNMerge'
	),
	function(
		x,
		params,
		...
	){

		stopifnot(valid(x, params))

		features <- SelectIntegrationFeatures(
			object.list = x,
			nfeatures = params@nfeatures,
			verbose = FALSE
		)

		object.list <- lapply(x, as.SingleCellExperiment)
		object.list <- lapply(object.list, function(obj) obj[features, ])
		out <- do.call(what = batchelor::fastMNN, args = c(object.list, list(d = params@ndims, k = params@n_neighbors)))

		# though fastMNN returns a reconstructed expression matrix, we don't store them in order to reduce memory usage
		x <- Reduce('merge', x)
		x[[params@reduction_name]] <- CreateDimReducObject(
			embeddings = reducedDim(out, 'corrected'),
			key = params@reduction_key,
			assay = params@normalize@assay_name
		)
		x@assays[[params@normalize@assay_name]]@var.features <- features
		x
	}
)


setClass(
	'LigerMerge', 
	representation(
    nrep = "integer",
    lambda = "numeric"
	),
	contains = c('BaseMerge'),
  prototype(
		nrep = 3L,
		lambda = 5,
		name = 'LigerMerge',
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0'),
			new('RPackage', package_name = 'rliger', package_version = '1.0.0')
		)
	)
)

#' Merge Seurat objects
#'
#' @param x SeuratList object containing single cell counts matrix
#' @param params a LigerMerge object 
#' @param ... Additional arguments
#' @return returns a Seurat object with integrated data
#' @importFrom Seurat GetAssayData
#' @export
#'
setMethod(
	'Merge',
	signature(
		x = 'SeuratList',
		params = 'LigerMerge'
	),
	function(
		x,
		params,
		...
	){

		stopifnot(valid(x, params))

		counts <- lapply(1:length(x), function(i){
			GetAssayData(
	 	    object = x[[i]],
				slot = 'counts',
				assay = params@normalize@assay_name 
			)
		})
		names(counts) <- names(x)

		features <- SelectIntegrationFeatures(
			object.list = x,
			nfeatures = params@nfeatures,
			verbose = FALSE
		)

		# Here we used the rliger pipeline for noramlization and scaling
		# Using scaled data from the univeral pipeline (e.g. Seurat) will sometimes have error 
		# see https://github.com/welch-lab/liger/issues/174
		# However, we still use the HVGs selected by the univeral pipeline
		#
		d <- rliger::createLiger(counts, remove.missing = FALSE)
		d <- rliger::normalize(d, remove.missing = FALSE)
		d@var.genes <- features
		d <- rliger::scaleNotCenter(d)
		d <- rliger::optimizeALS(d, k = params@ndims, verbose = FALSE)
		d <- rliger::quantile_norm(d, verbose = FALSE)

		x <- Reduce('merge', x)

		x[[params@reduction_name]] <- CreateDimReducObject(
			embeddings = d@H.norm,
			assay =  params@normalize@assay_name,
			key = params@reduction_key
		)
		x@assays[[params@normalize@assay_name]]@var.features <- features
		x
	}
)


setClass(
	'BBKNNMerge', 
	representation(
	),
	contains = 'BaseMerge', 
  prototype(
		name = 'BBKNNMerge',
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0'),
			new('PythonPackage', package_name = 'anndata', package_version = '0.7.8'),
			new('PythonPackage', package_name = 'scanpy', package_version = '1.8'),
			new('PythonPackage', package_name = 'bbknn', package_version = '1.5.1'),
			new('PythonPackage', package_name = 'leidenalg', package_version = '0.8.8'),
			new('RPackage', package_name = 'zellkonverter', package_version = '1.4.0'),
			new('RPackage', package_name = 'basilisk', package_version = '1.6.0'),
			new('PythonPackage', package_name = 'umap-learn', package_version = '0.5.2')
		)
	)
)

#' Merge Seurat objects by BBKNN (https://scanpy-tutorials.readthedocs.io/en/latest/integrating-data-using-ingest.html#BBKNN)
#'
#' @param x SeuratList object containing single cell counts matrix
#' @param params a BBKNNMerge object 
#' @param ... Additional arguments
#' @return returns a Seurat object with integrated data
#' @importFrom reticulate import
#' @importFrom Seurat RunPCA  CreateDimReducObject DefaultAssay GetAssayData
#' @export
#'
setMethod(
	'Merge',
	signature(
		x = 'SeuratList',
		params = 'BBKNNMerge'
	),
	function(
		x,
		params,
		...
	){

		stopifnot(valid(x, params))

	  bbknn <- import("bbknn")
		sc <- import("scanpy")
		anndata <- import("anndata")

		features <- SelectIntegrationFeatures(
			object.list = x,
			nfeatures = params@nfeatures,
			verbose = FALSE
		)

		x <- Reduce('merge', x)

		adata <- anndata$AnnData(
			X = t(GetAssayData(x, 'data')), 
			obs = x[[params@normalize@preprocess@batch]]
		)

		sc$tl$pca(adata, n_comps = params@ndims)

		sc$external$pp$bbknn(
			adata, 
			batch_key = params@normalize@preprocess@batch
		)

		sc$tl$umap(adata)
		y <- adata$obsm[['X_umap']]
		rownames(y) <- colnames(x)

		x[[params@reduction_name]] <- CreateDimReducObject(
			embeddings = y,
			assay =  params@normalize@assay_name,
			key = params@reduction_key
		)
		x@assays[[params@normalize@assay_name]]@var.features <- features
		x
	}
)


setClass(
	'ScanoramaMerge',
	representation(
	),
	contains = c('BaseMerge'),
  prototype(
		name = "ScanoramaMerge",
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0'),
			new('PythonPackage', package_name = 'anndata', package_version = '0.7.8'),
			new('PythonPackage', package_name = 'scanpy', package_version = '1.8'),
			new('PythonPackage', package_name = 'scanorama', package_version = '1.7.1'),
			new('RPackage', package_name = 'zellkonverter', package_version = '1.4.0'),
			new('RPackage', package_name = 'basilisk', package_version = '1.6.0')
		)
	)
)

#' Merge Seurat objects
#'
#' @param x SeuratList object containing single cell counts matrix
#' @param params a ScanoramaMerge object 
#' @param ... Additional arguments
#' @return returns a Seurat object with integrated data
#' @importFrom reticulate import py_capture_output
#' @importFrom Seurat SelectIntegrationFeatures VariableFeatures SplitObject GetAssayData CreateDimReducObject DefaultAssay
#' @export
#'
setMethod(
	'Merge',
	signature(
		x = 'SeuratList',
		params = 'ScanoramaMerge'
	),
	function(
		x,
		params,
		...
	){

		stopifnot(valid(x, params))
	  sc <- import("scanpy")
 	 	sr <- import("scanorama")
 	 	anndata <- import("anndata")

		for (i in 1:length(x)){
			x[[i]]@active.assay <- params@normalize@assay_name
		}

		features <- SelectIntegrationFeatures(
			object.list = x,
			nfeatures = params@nfeatures,
			verbose = FALSE
		)

		object.list <- lapply(x, function(obj) obj[features, ])

		assay_list <- list()
		gene_list <- list()
		for (i in 1:length(object.list)){
			assay_list[[i]] <- as.matrix(t(GetAssayData(object.list[[i]], "data")))
			gene_list[[i]] <- rownames(object.list[[i]])
		}

		integrated.corrected.data <- sr$correct(
			assay_list, 
			gene_list, 
			return_dimred = TRUE, 
			return_dense = TRUE, 
			dimred = params@ndims
		)

		intdimred <- do.call('rbind', integrated.corrected.data[[1]])
		rownames(intdimred) <- unlist(sapply(object.list, colnames))

		x <- Reduce('merge', x)

		x[[params@reduction_name]] <- CreateDimReducObject(
			embeddings = intdimred,
			key = params@reduction_key,
			assay = params@normalize@assay_name
		)
		x@assays[[params@normalize@assay_name]]@var.features <- features
		x	
	}
)


setClass(
	'scVIMerge',
	representation(
		n_layers = 'integer',								 
		encode_covariates = 'logical',
		deeply_inject_covariates = 'logical',
		use_layer_norm = "character",
		use_batch_norm = "character",
		max_epochs = 'integer',
		early_stopping = 'list',
		use_gpu = 'logical'
	),
	contains = c('BaseMerge'),
  prototype(
		name = "scVI",
		n_layers = 2L,
		encode_covariates = TRUE,
		deeply_inject_covariates = FALSE,
		use_layer_norm = "both",
		use_batch_norm = "none",
		max_epochs = 500L,
		early_stopping = list(
			early_stopping_metric = 'elbo',
			"save_best_state_metric"= "elbo",
			"patience"= 10,
			"threshold"= 0,
			"reduce_lr_on_plateau"= TRUE,
			"lr_patience"= 8,
			"lr_factor"= 0.1
		),
		use_gpu = TRUE,
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0'),
			new('PythonPackage', package_name = 'anndata', package_version = '0.7.8'),
			new('PythonPackage', package_name = 'scanpy', package_version = '1.8'),
			new('RPackage', package_name = 'zellkonverter', package_version = '1.4.0'),
			new('RPackage', package_name = 'basilisk', package_version = '1.6.0'),
			new('PythonPackage', package_name = 'scvi-tools', package_version = '0.14.5')
		)
	)
)

#' Merge Seurat objects by using scVI
#'
#' @param x SeuratList object containing single cell counts matrix
#' @param params a scVIMerge object 
#' @param ... Additional arguments
#' @return returns a Seurat object with integrated data
#' @importFrom Seurat VariableFeatures CreateDimReducObject DefaultAssay
#' @importFrom reticulate import
#' @export
#' @references https://scarches.readthedocs.io/en/latest/reference_building_from_scratch.html
#' @references Lopez, R., Regier, J., Cole, M.B. et al. Deep generative modeling for single-cell transcriptomics. Nat Methods 15, 1053–1058 (2018). https://doi.org/10.1038/s41592-018-0229-2
#'
setMethod(
	'Merge',
	signature(
		x= 'SeuratList',
		params = 'scVIMerge'
	),
	function(
		x,
		params,
		...
	){

		stopifnot(valid(x, params))

	  scvi <- import('scvi')
	 	anndata <- import("anndata")

		features <- SelectIntegrationFeatures(
			object.list = x,
			nfeatures = params@nfeatures,
			verbose = FALSE
		)

		x <- Reduce('merge', x)

		adata <- anndata$AnnData(
			X = t(GetAssayData(x, 'counts')[features, ]), 
			obs = x[[params@normalize@preprocess@batch]]
		)

		model <- .train_scvi_model(adata, params)

		model$get_normalized_expression()
		corrected <- model$get_normalized_expression()
		colnames(corrected) <- features
		x[[params@assay_name]] <- CreateAssayObject(
 		data = t(as(as.matrix(corrected), 'sparseMatrix'))
		)	

		latent <- model$get_latent_representation()
		rownames(latent) <- colnames(x)
		x[[params@reduction_name]] <- CreateDimReducObject(
			embeddings = latent,
			assay =  params@normalize@assay_name,
			key = params@reduction_key
		)
		x@assays[[params@normalize@assay_name]]@var.features <- features
		
		raw <-  GetAssayData(object = x, slot = 'counts')

                x@assays[[params@assay_name]]@counts <- raw[features,]

		x	
	}
)

.train_scvi_model <- function(adata, params){
	scvi <- import('scvi')
  scvi$settings$seed <- params@seed
 	scvi$model$SCVI$setup_anndata(adata, batch_key = params@normalize@preprocess@batch)
	model <- scvi$model$SCVI(
		adata, 
		n_latent = params@ndims,
		n_layers = params@n_layers,
		encode_covariates = params@encode_covariates,
		deeply_inject_covariates = params@deeply_inject_covariates,
		use_layer_norm = params@use_layer_norm,
		use_batch_norm = params@use_batch_norm
	)
	model$train(
		max_epochs = params@max_epochs,
		early_stopping = params@early_stopping,
		use_gpu = params@use_gpu
	)
	model
}


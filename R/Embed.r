#' Embed a SeuratList
#'
#' @param x a SeuratList object
#' @param params a PCAEmbed object
#' @param ... Additional arguments
#' @return returns a SeuratList object with latent embedding
#'
setMethod(
	'Embed',
	signature(
		x = 'SeuratList',
		params = 'BaseEmbed'
	),
	function(
		x,
		params,
		...
	){

		for (i in 1:length(x)){
			x[[i]] <- Embed(x[[i]], params, ...)
		}
		x
	}
)

#' The PCAEmbed class
#'
setClass(
	'PCAEmbed',
	representation(
	),
	contains = c('BaseEmbed'),
	prototype(
		name = 'PCAEmbed'
	)
)


#' Embed a Seurat object with PCA
#'
#' @param x a Seurat object
#' @param params a PCAEmbed object
#' @param ... Additional arguments
#' @return returns a data object with PCA embedding
#' @importFrom methods is
#' @importFrom Seurat ScaleData RunPCA
#'
setMethod(
	'Embed',
	signature(
		x = 'Seurat',
		params = 'PCAEmbed'
	),
	function(
		x,
		params,
		...
	){
		# to be implemented
		# 1. params@ndims should not be greater than # cells
		# 2. check whether the data has been preprocessed
		# stopifnot(valid(x, params))	

		x <- RunPCA(
			x,
			npcs = params@ndims,
			reduction.name = params@reduction_name,
			reduction.key = params@reduction_key,
			verbose = FALSE
		)
		x
	}
)


#' The scCCESSSEmbed class
#'
setClass(
	'scCCESSSEmbed',
	representation(
		batch_size = 'integer',
		max_random_projection = 'integer',
		hidden_dims = 'integer',
		learning_rate = 'numeric',
		epochs = 'integer'
	),
	contains = c('BaseEmbed'),
	prototype(
		name = 'scCCESSSEmbed',
		batch_size = 64L,
		max_random_projection = 2048L,
		hidden_dims = 128L,
		learning_rate = 0.001,
		epochs = 100L,
		dependences = list(
			new('RPackage', package_name = 'tensorflow', package_version = '2.8.0'),
			new('RPackage', package_name = 'keras', package_version = '2.8.0'),
			new('RPackage', package_name = 'clue', package_version = '0.3.60'),
			new('RPackage', package_name = 'scCCESS', package_version = '0.2.0'),
			new('PythonPackage', package_name = 'tensorflow', package_version = '2.8.0'),
			new('PythonPackage', package_name = 'keras', package_version = '2.8.0')
		)
	)
)


#' Embed a Seurat object with scCCESSS
#'
#' @param x a Seurat object
#' @param params a scCCESSSEmbed object
#' @param ... Additional arguments
#' @return returns a data object with embedding results
#' @importFrom Seurat CreateDimReducObject 
#'
setMethod(
	'Embed',
	signature(
		x = 'Seurat',
		params = 'scCCESSSEmbed'
	),
	function(
		x,
		params,
		...
	){
		# to be implemented
		# 1. params@ndims should not be greater than # cells
		# 2. check whether the data has been preprocessed
		# stopifnot(valid(x, params))	

		raw_assay <- params@preprocess@raw_assay

		results <- scCCESS::encode(
			x@assays[[raw_assay]]@scale.data,
			seed = params@seed,
			max_random_projection = params@max_random_projection,
			encoded_dim = params@ndims,
			hidden_dims = params@hidden_dims,
			learning_rate = params@learning_rate,
			batch_size = params@batch_size,
			epochs = params@epochs,
			verbose = 0, 	# silent
			scale = FALSE, 
			genes_as_rows = TRUE
		)

		x[[params@name]] <- CreateDimReducObject(
			embeddings = results,
			key = params@reduction_key,
			assay = DefaultAssay(x)
		)

		x
	}
)


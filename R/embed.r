#' Embed a SeuratList
#'
#' @param x a SeuratList object
#' @param params a PCAEmbed object
#' @param ... Additional arguments
#' @export
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

setClass(
	'PCAEmbed',
	representation(
	),
	contains = c('BaseEmbed'),
	prototype(
		name = 'PCAEmbed'
	)
)

setMethod('initialize', 'PCAEmbed', function(.Object, check_dependencies = TRUE, ...){
	callNextMethod(.Object, check_dependencies = check_dependencies, ...)
})


#' Embed a Seurat object with PCA
#'
#' @param x a Seurat object
#' @param params a PCAEmbed object
#' @param ... Additional arguments
#' @return returns a data object with PCA embedding
#' @importFrom methods is
#' @importFrom Seurat ScaleData RunPCA
#' @export
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
		x <- RunPCA(
			x,
			assay = params@normalize@assay_name,
			npcs = params@ndims,
			rev.pca = FALSE,
			weight.by.var = TRUE,
			reduction.name = params@reduction_name,
			reduction.key = params@reduction_key,
			seed.use = params@seed,
			approx = TRUE,
			verbose = FALSE
		)
		x
	}
)


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
#' @export
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

		raw_assay <- params@normalize@assay_name

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


setClass(
	'UMAPEmbed',
	representation(
		embedding = 'BaseEmbed',
		n_neighbors = 'integer',
		metric = 'character'
	),
	contains = c('BaseEmbed'),
	prototype(
		n_neighbors = 30L,
		ndims = 2L,
		metric = 'cosine'
	)
)

#' @importFrom methods callNextMethod
#'
setMethod('initialize', 'UMAPEmbed', function(.Object, check_dependencies = TRUE, ...){
	.Object <- callNextMethod(.Object, check_dependencies = check_dependencies, ...)
	.Object@name <- sprintf('%sUMAP', .Object@embedding@name)
	callNextMethod(.Object, check_dependencies = check_dependencies, ...)
})


#' Embed a Seurat object with UMAP
#'
#' @param x a Seurat object
#' @param params a UMAPEmbed object
#' @param ... Additional arguments
#' @return returns a data object with embedding results
#' @importFrom Seurat CreateDimReducObject RunUMAP
#' @export
#'
setMethod(
	'Embed',
	signature(
		x = 'Seurat',
		params = 'UMAPEmbed'
	),
	function(
		x,
		params,
		...
	){
		x <- RunUMAP(
			x,
			dims = 1:params@embedding@ndims,
			reduction = params@embedding@reduction_name,
			assay = params@normalize@assay_name,
			slot = 'data',
			umap.method = "uwot",
			n.neighbors = params@n_neighbors,
			n.components = params@ndims,
			metric = params@metric,
			seed.use = params@seed,
			verbose = FALSE,
			reduction.name = params@reduction_name,
			reduction.key = params@reduction_key
		)
		x
	}
)


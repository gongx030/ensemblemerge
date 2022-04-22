setClass(
	'symphonyReferenceMap',
	representation(
		k = 'integer',
		vargenes_method = 'character',
		theta = 'numeric',
		ndims = 'integer',
		sigma = 'numeric'
	),
	contains = c('BaseReferenceMap'),
	prototype(
		name = 'symphonyReferenceMap',
		k = 100L,
		vargenes_method = 'vst',
		theta = 2,
		ndims = 20L,
		sigma = 0.1,
		dependences = list(
			new('RPackage', package_name = 'symphony', package_version = '0.1.0')
		)
	)
)


setMethod('initialize', 'symphonyReferenceMap', function(.Object, ...){
	callNextMethod(.Object, ...)
})


#' Map a SeuratList object (query) onto another SeruatList object (atlas) by Symphony (https://github.com/immunogenomics/symphony)
#'
#' @param query a SeuratList object with multiple batches
#' @param atlas a SeuratList object with multiple batches
#' @param params a symphonyReferenceMap object
#' @param ... Additional arguments
#' @return returns a SeuratList object with query and atlas data with integrated dimension reduction
#' @export
#'
setMethod(
	'ReferenceMap',
	signature(
		query = 'SeuratList',
		atlas = 'SeuratList',
		params = 'symphonyReferenceMap'
	),
	function(
		query,
		atlas,
		params,
		...
	){
		query <- Reduce('merge', query)
		atlas <- Reduce('merge', atlas)
		ReferenceMap(query, atlas, params)
	}
)


#' Map a Seurat object (query) onto a SeruatList object (atlas) by Symphony (https://github.com/immunogenomics/symphony)
#'
#' @param query a Seurat object
#' @param atlas a SeuratList object with multiple batches
#' @param params a symphonyReferenceMap object
#' @param ... Additional arguments
#' @return returns a SeuratList object with query and atlas data with integrated dimension reduction
#' @export
#'
setMethod(
	'ReferenceMap',
	signature(
		query = 'Seurat',
		atlas = 'SeuratList',
		params = 'symphonyReferenceMap'
	),
	function(
		query,
		atlas,
		params,
		...
	){
		atlas <- Reduce('merge', atlas)
		ReferenceMap(query, atlas, params)
	}
)


#' Map a SeuratList object (query) onto another SeruatList object (atlas) by Symphony (https://github.com/immunogenomics/symphony)
#'
#' @param query a SeuratList object with multiple batches
#' @param atlas a Seurat object
#' @param params a symphonyReferenceMap object
#' @param ... Additional arguments
#' @return returns a SeuratList object with query and atlas data with integrated dimension reduction
#' @export
#'
setMethod(
	'ReferenceMap',
	signature(
		query = 'SeuratList',
		atlas = 'Seurat',
		params = 'symphonyReferenceMap'
	),
	function(
		query,
		atlas,
		params,
		...
	){
		query <- Reduce('merge', query)
		ReferenceMap(query, atlas, params)
	}
)


#' Map a Seurat object (query) onto another Seruat object (atlas) by Symphony (https://github.com/immunogenomics/symphony)
#'
#' @param query a Seruat object
#' @param atlas a Seurat object
#' @param params a symphonyReferenceMap object
#' @param ... Additional arguments
#' @return returns a SeuratList object with query and atlas data with integrated dimension reduction
#' @references Kang, J.B., Nathan, A., Weinand, K. et al. Efficient and precise single-cell reference atlas mapping with Symphony. Nat Commun 12, 5890 (2021). https://doi.org/10.1038/s41467-021-25957-x
#' @export
#'
setMethod(
	'ReferenceMap',
	signature(
		query = 'Seurat',
		atlas = 'Seurat',
		params = 'symphonyReferenceMap'
	),
	function(
		query,
		atlas,
		params,
		...
	){

		if (any(colnames(query) %in% colnames(atlas)))
			stop('there are duplicated cell names between query and atlas')

		atlas@active.assay <- params@normalize_atlas@assay_name

		vars_atlas <- NULL
		if (length(params@normalize_atlas@preprocess@batch) > 0)
			vars_atlas <-  params@normalize_atlas@preprocess@batch

		# Build reference
		reference <- symphony::buildReference(
			GetAssayData(atlas, 'data'),                   # reference expression (genes by cells)
			atlas@meta.data,              # reference metadata (cells x attributes)
			vars = vars_atlas,         # variable(s) to integrate over
			K = params@k,                   # number of Harmony soft clusters
			verbose = TRUE,            # display verbose output
			do_umap = FALSE,            # run UMAP and save UMAP model to file
			do_normalize = FALSE,      # perform log(CP10k) normalization on reference expression
			vargenes_method = params@vargenes_method,   # variable gene selection method: 'vst' or 'mvp'
			vargenes_groups = NULL, # metadata column specifying groups for variable gene selection within each group
			topn = params@normalize_atlas@numHVG,               # number of variable genes (per group)
			theta = params@theta,                 # Harmony parameter(s) for diversity term
			d = params@ndims,                    # number of dimensions for PCA
			save_uwot_path = NULL, # file path to save uwot UMAP model
			additional_genes = NULL    # vector of any additional genes to force include
		)

		query@active.assay <- params@normalize_query@assay_name

		vars_query <- NULL
		if (length(params@normalize_query@preprocess@batch) > 0)
			vars_query <-  params@normalize_query@preprocess@batch

		results <- symphony::mapQuery(
			GetAssayData(query, 'data'),             # query gene expression (genes x cells)
			query@meta.data,        # query metadata (cells x attributes)
			reference,             # Symphony reference object
			vars = vars_query,           # Query batch variables to harmonize over (NULL treats query as one batch)
			do_normalize = FALSE,  # perform log(CP10k) normalization on query (set to FALSE if already normalized)
			do_umap = FALSE,
			sigma = 0.1
		)  

		query[[params@reduction_name]] <- CreateDimReducObject(
			embeddings = results$Z %>% t(),
			key = params@reduction_key,
			assay = params@normalize_query@assay_name
		)

		atlas[[params@reduction_name]] <- CreateDimReducObject(
			embeddings = reference$Z_corr %>% t(),
			key = params@reduction_key,
			assay = params@normalize_atlas@assay_name
		)

		new('SeuratList', list(query = query, atlas = atlas))
	}
)


setClass(
	'SeuratReferenceMap',
	representation(
		normalization.method = 'character',
		reduction_type = 'character',
		reference_reduction = 'character'
	),
	contains = c('BaseReferenceMap'),
	prototype(
		name = 'SeuratReferenceMap',
		ndims = 20L,
		reduction_type = 'pcaproject',
		dependences = list(
			new('RPackage', package_name = 'Seurat', package_version = '4.1.0')
		)
	),
	validity = function(object){
		msg <- NULL
		if (class(object@normalize_query) != class(object@normalize_atlas))
			msg <- 'The class of object@normalize_query and object@normalize_atlas must be the same'
		return(msg)
	}
)

setMethod('initialize', 'SeuratReferenceMap', function(.Object, ...){
	.Object <- callNextMethod(.Object, ...)
	if (is(.Object@normalize_atlas, 'SeuratNormalize'))
		.Object@normalization.method <- 'LogNormalize'
	else if (is(.Object@normalize_atlas, 'SCTransformNormalize'))
		.Object@normalization.method <- 'SCT'
	.Object
})


#' Map a SeuratList object (query) onto another Seurat object (atlas) by Seurat (https://satijalab.org/seurat/articles/integration_mapping.html)
#'
#' @param query a SeuratList object with multiple batches
#' @param atlas a SeuratList object with multiple batches
#' @param params a SeuratReferenceMap object
#' @param ... Additional arguments
#' @return returns a SeuratList object with query and atlas data with integrated dimension reduction
#' @export
#'
setMethod(
	'ReferenceMap',
	signature(
		query = 'SeuratList',
		atlas = 'SeuratList',
		params = 'SeuratReferenceMap'
	),
	function(
		query,
		atlas,
		params,
		...
	){

		params_merge <- new('SeuratMerge', normalize = params@normalize_query, ndims = params@ndims)
		query <- Merge(query, params_merge)

		params_merge <- new('SeuratMerge', normalize = params@normalize_atlas, ndims = params@ndims)
		atlas <- Merge(atlas, params_merge)

		params@reference_reduction <- params_merge@reduction_name

		.seurat_referencemap_core(query, atlas, params, ...)

	}
)



#' Map a Seurat object (query) onto another Seurat object (atlas) by Seurat (https://satijalab.org/seurat/articles/integration_mapping.html)
#'
#' @param query a Seurat object 
#' @param atlas a Seurat object 
#' @param params a SeuratReferenceMap object
#' @param ... Additional arguments
#' @return returns a SeuratList object with query and atlas data with integrated dimension reduction
#' @export
#'
setMethod(
	'ReferenceMap',
	signature(
		query = 'Seurat',
		atlas = 'Seurat',
		params = 'SeuratReferenceMap'
	),
	function(
		query,
		atlas,
		params,
		...
	){

		params_pca <- new('PCAEmbed', normalize = params@normalize_atlas, ndims = params@ndims)
		atlas <- Embed(atlas, params_pca)

		params@reference_reduction <- params_pca@reduction_name

		.seurat_referencemap_core(query, atlas, params, ...)
	}
)


#' Map a SeuratList object (query) onto another Seurat object (atlas) by Seurat (https://satijalab.org/seurat/articles/integration_mapping.html)
#'
#' @param query a SeuratList object (with batches)
#' @param atlas a Seurat object 
#' @param params a SeuratReferenceMap object
#' @param ... Additional arguments
#' @return returns a SeuratList object with query and atlas data with integrated dimension reduction
#' @export
#'
setMethod(
	'ReferenceMap',
	signature(
		query = 'SeuratList',
		atlas = 'Seurat',
		params = 'SeuratReferenceMap'
	),
	function(
		query,
		atlas,
		params,
		...
	){

		params_merge <- new('SeuratMerge', normalize = params@normalize_query, ndims = params@ndims)
	  query <- Merge(query, params_merge)

		params_pca <- new('PCAEmbed', normalize = params@normalize_atlas, ndims = params@ndims)
		atlas <- Embed(atlas, params_pca)
		params@reference_reduction <- params_pca@reduction_name

		.seurat_referencemap_core(query, atlas, params, ...)
	}
)

#' Map a Seurat object (query) onto another SeuratList object (atlas) by Seurat (https://satijalab.org/seurat/articles/integration_mapping.html)
#'
#' @param query a Seurat object
#' @param atlas a SeuratList object  (with batches)
#' @param params a SeuratReferenceMap object
#' @param ... Additional arguments
#' @return returns a SeuratList object with query and atlas data with integrated dimension reduction
#' @export
#'
setMethod(
	'ReferenceMap',
	signature(
		query = 'Seurat',
		atlas = 'SeuratList',
		params = 'SeuratReferenceMap'
	),
	function(
		query,
		atlas,
		params,
		...
	){

		params_merge <- new('SeuratMerge', normalize = params@normalize_atlas, ndims = params@ndims)
		atlas <- Merge(atlas, params_merge)
		params@reference_reduction <- params_merge@reduction_name

		.seurat_referencemap_core(query, atlas, params, ...)
	}
)


.seurat_referencemap_core <- function(query, atlas, params, ...){

	x <- FindTransferAnchors(
		reference = atlas, 
		query = query,
		normalization.method = params@normalization.method,
		reference.assay = NULL,
		reference.neighbors = NULL,
		query.assay = NULL,
		reduction = params@reduction_type,
		reference.reduction = params@reference_reduction,
		project.query = FALSE,
		features = NULL,
		scale = FALSE,
		npcs = params@ndims,
		l2.norm = TRUE,
		dims = 1:params@ndims,
		k.anchor = 5,
		k.filter = 200,
		k.score = 30,
		max.features = 200,
		nn.method = "annoy",
		n.trees = 50,
		eps = 0,
		approx.pca = TRUE,
		mapping.score.k = NULL,
		verbose = TRUE
	)

	x <- IntegrateEmbeddings(
		x,
		reference = atlas,
		query = query,
		new.reduction.name = params@reduction_name,
		reductions = params@reduction_type,
		dims.to.integrate = NULL,
		k.weight = 100,
		weight.reduction = NULL,
		reuse.weights.matrix = TRUE,
		sd.weight = 1,
		preserve.order = FALSE,
		verbose = TRUE
	)

	atlas[[params@reduction_name]] <- CreateDimReducObject(
		embeddings = atlas[[params@reference_reduction]]@cell.embeddings,
		key = params@reduction_key,
		assay = params@normalize_atlas@assay_name
	)

	new('SeuratList', list(query = x, atlas = atlas))
}

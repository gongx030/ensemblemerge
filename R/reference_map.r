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


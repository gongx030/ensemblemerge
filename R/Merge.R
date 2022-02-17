#' Merge SummarizedExperiment objects
#'
#' @param params a BaseMerge object 
#' @param data a SummarizedExperiment object 
#' @param ... Additional arguments
#'
#' @return returns a data object with integrated data
#'
#' @importFrom Seurat as.SingleCellExperiment
#'
setMethod(
	'Merge',
	signature(
		params = 'BaseMerge',
		data  = 'SummarizedExperiment'
	),
	function(
		params,
		data,
		...
	){

		data <- as.SingleCellExperiment(data)
		Merge(params, data, ...)

})


#' Merge Seurat objects
#'
#' @param params a BaseMerge object 
#' @param data Seurat object containing single cell counts matrix
#' @param ... Additional arguments
#' @return returns a data object with integrated data
#' @importFrom Seurat as.SingleCellExperiment
#' @importFrom methods is
#'
setMethod(
	'Merge',
	signature(
		params = 'BaseMerge',
		data  = 'Seurat'
	),
	function(
		params,
		data,
		...
	){

		if (is(params, 'SeuratParams')){
    	data = run_Seurat(params, data)
		}else if (is(params, 'HarmonyParams')){
    	data = run_Harmony(params, data)
		}else if (is(params, 'FastMNNParams')){
    	data = run_fastMNN(params, data)
		}else if (is(params, 'UncorrectedParams')){
	    data = run_Uncorrected(params, data)
		}else if (is(params, 'LigerParams')){
	    data = run_Liger(params, data)
		}else if (is(params, 'BBKNNParams')){
	    data = run_BBKNN(params, data)
		}else if (is(params, 'ScanoramaParams')){
	    data = run_Scanorama(params, data)
		}else if (is(params, 'scVIParams')){
	    data = run_scVI(params, data)
		}else
			stop(sprintf('unknown params class: %s', class(params)))

		data
})

#' Merge SummarizedExperiment objects
#'
#' @param params a ParamsList object 
#' @param data a SummarizedExperiment object 
#' @param ... Additional arguments
#'
#' @importFrom Seurat CreateSeuratObject 
#' @importFrom SummarizedExperiment assays colData colData<- rowData rowData<- 
#'
#' @return returns a data object with integrated data
#'
#' @export
#'
setMethod(
	'Merge',
	signature(
		params = 'ParamsList',
		data  = 'SummarizedExperiment'
	),
	function(
		params,
		data,
		counts = 'counts',
		assay = 'RNA',
		...
	){

		seurat <- CreateSeuratObject(
			counts = assays(data)[[counts]], # Unnormalized data such as raw counts or TPMs
			meta.data = as.data.frame(colData(data)),
			assay = assay
		)

		if (ncol(rowData(data)) > 0){
			seurat[[assaye]][[]] <- as.data.frame(rowData(data))
		}

		Merge(params, seurat, ...)
})

#' Merge Seurat objects
#'
#' @param params a ParamsList object 
#' @param data a Seurat object 
#' @param ... Additional arguments
#'
#' @return returns a data object with integrated data
#'
#' @export
#'
setMethod(
	'Merge',
	signature(
		params = 'ParamsList',
		data  = 'Seurat'
	),
	function(
		params,
		data,
		...
	){

		stopifnot(!any(duplicated(colnames(data))))

		# use the preprocessing pipeline from the first one to process the data
		# need to deal with the situation where the preprocessing specificication may be different from constituent methods
		# we should use consistent preprocessing method for constituent method
		#
		data <- Preprocess(params[[1L]], data)

		for (i in 1:length(params)){
			data <- Merge(params[[i]], data)
		}
		data
})


#' Merge Seurat objects by the EnsembleMerge method
#'
#' @param params a EnsembleMergeParams object
#' @param data a Seurat object 
#' @param ... Additional arguments
#' @return returns a data object (SingleCellExperiment or Seurat) of the integrated data
#' @importFrom uwot umap
#' @importFrom Seurat DefaultAssay
#' @export
#'
setMethod(
	'Merge',
	signature(
		params = 'EnsembleMergeParams',
		data  = 'Seurat'
	),
	function(
		params,
		data,
		...
	){

		integrated <- list()

		for (i in 1:length(params@constituent)){
			data <- Merge(params@constituent[[i]], data)
		}

		browser()

		for (i in 1:length(params@constituent)){
			integrated[[i]] <- Merge(params@constituent[[i]], data)
		}

		cn <- colnames(data)
		res <- ensemblemerge_core(params, integrated)
		ng <- res$ng[cn, cn]
		y <- umap(ng)
		rownames(y) <- colnames(data)
		data[[params@umap_name]] <- CreateDimReducObject(
			embeddings = y, 
			key = params@umap_key, 
			assay = DefaultAssay(data)
		)
		data@misc$kta_weight <- res$weight
		data
	}
)

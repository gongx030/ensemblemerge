#' Merge SummarizedExperiment objects
#'
#' @param data a SummarizedExperiment object 
#' @param params a BaseMerge object 
#' @param ... Additional arguments
#'
#' @return returns a data object with integrated data
#'
#' @importFrom Seurat as.SingleCellExperiment
#'
setMethod(
	'Merge',
	signature(
		data  = 'SummarizedExperiment',
		params = 'BaseMerge'
	),
	function(
		data,
		params,
		...
	){

		data <- as.SingleCellExperiment(data)
		Merge(params, data, ...)

})


#' Merge Seurat objects
#'
#' @param data Seurat object containing single cell counts matrix
#' @param params a BaseMerge object 
#' @param ... Additional arguments
#' @return returns a data object with integrated data
#' @importFrom Seurat as.SingleCellExperiment
#' @importFrom methods is
#'
setMethod(
	'Merge',
	signature(
		data  = 'Seurat',
		params = 'BaseMerge'
	),
	function(
		data,
		params,
		...
	){

		sprintf('Merge | running %s', params@name) %>% message()

		stopifnot(valid(data, params))

		if (is(params, 'SeuratMerge')){
    	data = run_Seurat(params, data)
		}else if (is(params, 'HarmonyMerge')){
    	data = run_Harmony(params, data)
		}else if (is(params, 'FastMNNMerge')){
    	data = run_fastMNN(params, data)
		}else if (is(params, 'UncorrectedMerge')){
	    data = run_Uncorrected(params, data)
		}else if (is(params, 'LigerMerge')){
	    data = run_Liger(params, data)
		}else if (is(params, 'BBKNNMerge')){
	    data = run_BBKNN(params, data)
		}else if (is(params, 'ScanoramaMerge')){
	    data = run_Scanorama(params, data)
		}else if (is(params, 'scVIMerge')){
	    data = run_scVI(params, data)
		}else
			stop(sprintf('unknown params class: %s', class(params)))

		sprintf('Merge | running UMAP') %>% message()
	  data <- RunUMAP(
			data, 
			reduction = params@name, 
			dims = 1:params@npcs, 
			reduction.key = params@umap_key, 
			seed.use = 1, 
			reduction.name = params@umap_name, 
			verbose = FALSE
		)

		data
})

#' Merge SummarizedExperiment objects
#'
#' @param data a SummarizedExperiment object 
#' @param params a MethodList object 
#' @param counts the assay field for raw counts in a SingleCellExperiment object (default: 'counts')
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
		data  = 'SummarizedExperiment',
		params = 'Params'
	),
	function(
		data,
		params,
		counts = 'counts',
		...
	){

		seurat <- CreateSeuratObject(
			counts = assays(data)[[counts]], # Unnormalized data such as raw counts or TPMs
			meta.data = as.data.frame(colData(data)),
			assay = params@preprocess@raw_assay
		)

		if (ncol(rowData(data)) > 0){
			seurat[[params@preprocess@raw_assay]][[]] <- as.data.frame(rowData(data))
		}

		Merge(seurat, params, ...)
})

#' Merge Seurat objects
#'
#' @param data a Seurat object 
#' @param params a MethodList object 
#' @param ... Additional arguments
#'
#' @return returns a data object with integrated data
#'
#' @export
#'
setMethod(
	'Merge',
	signature(
		data  = 'Seurat',
		params = 'Params'
	),
	function(
		data,
		params,
		...
	){

		# use the preprocessing pipeline from the first one to process the data
		# need to deal with the situation where the preprocessing specificication may be different from constituent methods
		# we should use consistent preprocessing method for constituent method
		#
		data <- Preprocess(data, params@preprocess)

		for (i in 1:length(params@constituent)){
			data <- Merge(data, params@constituent[[i]])
		}
		data
})



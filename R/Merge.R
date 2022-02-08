#' Merge SummarizedExperiment objects
#'
#' @param params a BaseMerge object 
#' @param data a SummarizedExperiment object 
#' @param ... Additional arguments
#'
#' @return returns a data object with integrated data
#'
#' @export
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

#' Merge SingleCellExperiment objects
#'
#' @param params a BaseMerge object 
#' @param data a SingleCellExperiment
#' @param counts The assays field that contains the raw counts data
#' @param ... Additional arguments
#'
#' @return returns a data object with integrated data
#'
#' @importFrom Seurat CreateSeuratObject 
#' @importFrom SummarizedExperiment assays colData colData<- rowData rowData<- 
#'
#' @export
#'
setMethod(
	'Merge',
	signature(
		params = 'BaseMerge',
		data  = 'SingleCellExperiment'
	),
	function(
		params,
		data,
		counts = 'counts',
		...
	){

		seurat <- CreateSeuratObject(
			counts = assays(data)[[counts]],
			meta.data = as.data.frame(colData(data))
		)

		if (ncol(rowData(data)) > 0){
			seurat[["RNA"]][[]] <- as.data.frame(rowData(data))
		}

		Merge(params, seurat, ...)
	}	
)

#' Merge Seurat objects
#'
#' @param params a BaseMerge object 
#' @param data Seurat object containing single cell counts matrix
#' @param ... Additional arguments
#' @return returns a data object with integrated data
#' @importFrom Seurat as.SingleCellExperiment
#' @importFrom methods is
#' @export
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

		stopifnot(!any(duplicated(colnames(data))))

		data <- Preprocess(params, data)

		if (is(params, 'SeuratParams')){
    	integrated = run_Seurat(params, data)
		}else if (is(params, 'HarmonyParams')){
    	integrated = run_Harmony(params, data)
		}else if (is(params, 'FastMNNParams')){
    	integrated = run_fastMNN(params, data)
		}else if (is(params, 'UncorrectedParams')){
	    integrated = run_Uncorrected(params, data)
		}else if (is(params, 'LigerParams')){
	    integrated = run_Liger(params, data)
		}else if (is(params, 'BBKNNParams')){
	    integrated = run_BBKNN(params, data)
		}else if (is(params, 'ScanoramaParams')){
	    integrated = run_Scanorama(params, data)
		}else if (is(params, 'scVIParams')){
	    integrated = run_scVI(params, data)
		}else
			stop(sprintf('unknown params class: %s', class(params)))

		if(params@return == "Seurat"){
			return(integrated)
		}else if(params@return == "SingleCellExperiment"){
			integrated = as.SingleCellExperiment(integrated)
			return(integrated)
		}else{
			stop("Invalid return type, check params@return")
		}
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

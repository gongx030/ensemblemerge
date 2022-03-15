#' Print class BasePreprocess
#' 
#' @param object a BasePreprocess object
#'
#' @importFrom methods is
#'
setMethod(
	'show',
	signature(
		object = 'BasePreprocess'
	),
	function(
		object
	){
		sprintf('Preprocessing:') %>% message()
		sprintf('%20.20s: %21.d [%s]', 'min_cells', object@min_cells, 'mininum number of cells a gene is expressed') %>% message()
		sprintf('%20.20s: %21.d [%s]', 'min_genes', object@min_genes, 'mininum number of expressed genes a cell should have') %>% message()
		sprintf('%20.20s: %21.21s [%s]', 'norm_data', object@norm_data, 'whether or not the data need to be normalized') %>% message()
		sprintf('%20.20s: %21.21s [%s]', 'norm_method', object@norm_method, 'the normalization method') %>% message()
		sprintf('%20.20s: %21.d [%s]', 'scale_factor', object@scale_factor, 'Scale factor') %>% message()
		sprintf('%20.20s: %21.21s [%s]', 'scaling', object@scaling, 'whether or data need to scaled') %>% message()
		sprintf('%20.20s: %21.d [%s]', 'numHVG', object@numHVG, 'number of highly variable features') %>% message()
		sprintf('%20.20s: %21.21s [%s]', 'raw_assay', object@raw_assay, 'the raw counts field in the Seurat object') %>% message()
		sprintf('%20.20s: %21.21s [%s]', 'batch', object@batch, 'Batch field name') %>% message()
	}
)

#' Print class SeuratPreprocess 
#' 
#' @param object a SeuratPreprocess object
#'
#' @importFrom methods is
#'
setMethod(
	'show',
	signature(
		object = 'SeuratPreprocess'
	),
	function(
		object
	){
		callNextMethod()
		sprintf('%20.20s: %21.21s [%s]', 'selection.method', object@selection.method, 'Features selection method') %>% message()
		sprintf('%20.20s: %21.21s [%s]', 'batchwise', object@batchwise, 'whether or not select features for each batch') %>% message()
	}
)

#' Print class BaseMerge
#' 
#' @param object a SeuratPreprocess object
#'
#' @importFrom methods is
#'
setMethod(
	'show',
	signature(
		object = 'BaseMerge'
	),
	function(
		object
	){
		sprintf('%s:', object@name) %>% message()
		sprintf('%20.20s: %21.d [%s]', 'npcs', object@npcs, 'Size of latent dimensions') %>% message()
		sprintf('%20.20s: %21.d [%s]', 'umap_dim', object@umap_dim, 'Size of UMAP dimensions') %>% message()
		sprintf('%20.20s: %21.21s [%s]', 'umap_name', object@umap_name, 'UMAP reduction') %>% message()
	}
)


#' Print class SeuratMerge
#' 
#' @param object a SeuratMerge object
#' @export
#'
setMethod(
	'show',
	signature(
		object = 'SeuratMerge'
	),
	function(
		object
	){
		callNextMethod()
		sprintf('%20.10s: %21.3f [%s]', 'k.weight', object@k.weight, 'weight for neighbor function') %>% message()
	}
)


#' Print class HarmonyMerge 
#' 
#' @param object a HarmonyMerge object
#' @export
#'
setMethod(
	'show',
	signature(
		object = 'HarmonyMerge'
	),
	function(
		object
	){
		callNextMethod()
		sprintf('%20.20s: %21.3f [%s]', 'theta', object@theta, 'diversity clustering penalty parameter, larger values increase diversity') %>% message()
		sprintf('%20.20s: %21.d [%s]', 'max_iter_cluster', object@max_iter_cluster, 'maximum number of learning iterations per cluster') %>% message()
	}
)

#' Print class FastMNNMerge 
#' 
#' @param object a FastMNNMerge object
#' @export
#'
setMethod(
	'show',
	signature(
		object = 'FastMNNMerge'
	),
	function(
		object
	){
		callNextMethod()
		sprintf('%20.20s: %21.d [%s]', 'n_neighbors', object@n_neighbors, 'number of neighbors used in calculating neighboring graph') %>% message()
	}
)

#' Print class LigerMerge 
#' 
#' @param object a LigerMerge object
#' @export
#'
setMethod(
	'show',
	signature(
		object = 'LigerMerge'
	),
	function(
		object
	){
		callNextMethod()
		sprintf('%20.20s: %21.d [%s]', 'nrep', object@nrep, 'number of repeats') %>% message()
		sprintf('%20.20s: %21.3f [%s]', 'lambda', object@lambda, 'lambda') %>% message()
	}
)

#' Print class BBKNNMerge 
#' 
#' @param object a BBKNNMerge object
#' @export
#'
setMethod(
	'show',
	signature(
		object = 'BBKNNMerge'
	),
	function(
		object
	){
		callNextMethod()
	}
)



#' Print class Seurat
#' 
#' @param object a SeuratPreprocess object
#'
#' @importFrom methods is
#'
setMethod(
	'show',
	signature(
		object = 'Params'
	),
	function(
		object
	){
		show(object@preprocess)
		for (i in 1:length(object@constituent)){
			sprintf('%s:', object@constituent[[i]]@name) %>% message()
			show(object@constituent[[i]])
		}
	}
)

#' Print class EnsembleMerge
#' 
#' @param object a EnsembleMerge object
#' @export
#'
setMethod(
	'show',
	signature(
		object = 'EnsembleMerge'
	),
	function(
		object
	){
		sprintf('%s:', object@name) %>% message()
		sprintf('%20.20s: %21.d [%s]', 'umap_dim', object@umap_dim, 'Size of UMAP dimensions') %>% message()
		sprintf('%20.20s: %21.21s [%s]', 'umap_name', object@umap_name, 'UMAP reduction') %>% message()
		sprintf('%20.20s: %21.d [%s]', 'k.param', object@k.param, 'k for KNN') %>% message()
		sprintf('%20.20s: %s', 'reductions', paste(object@constituent_reduction_names, collapse = ',')) %>% message()
	}
)

#' Get a neighbor graph
#' 
#' @param data a Seurat object
#' @param params a BaseMerge object
#' @param ... Additional arguments
#'
#' @importFrom Seurat FindNeighbors RunUMAP 
#'
#' @return returns a Seurat object with neighboring information
#' @export
#'
setMethod(
	"getNeighborGraph", 
	signature(
		data = 'Seurat',
		params = "BaseMerge"
	),
	function(data, params, ...) {

		if (params@latent){
    	data <- FindNeighbors(
				data, 
				compute.SNN = TRUE, 
				reduction = params@name,
				dims = 1:params@npcs, 
				graph.name = c(params@knn_name, params@snn_name),
			)
		}else{
			data <- RunUMAP(
				data, 
				reduction = params@name, 
				dims = 1:params@npcs,
				reduction.name = params@umap_name, 
				reduction.key = params@umap_key
			)
    	data <- FindNeighbors(
				data, 
				compute.SNN = TRUE, 
				reduction = params@umap_name, 
				dims = 1:params@umap_dim, 
				graph.name = c(params@knn_name, params@snn_name)
			)	
		}
		data
	}
)

#' Get a neighbor graph
#' 
#' @param data a Seurat object
#' @param params a EnsembleMerge object
#' @param ... Additional arguments
#'
#' @importFrom Seurat FindNeighbors 
#'
#' @return returns a Seurat object with neighboring information
#' @export
#'
setMethod(
	"getNeighborGraph", 
	signature(
		data = 'Seurat',
		params = "EnsembleMerge"
	),
	function(data, params, ...) {

   	data <- FindNeighbors(
			data, 
			compute.SNN = TRUE, 
			reduction = params@umap_name,
			dims = 1:params@umap_dim, 
			graph.name = c(params@knn_name, params@snn_name)
		)
		data
	}
)

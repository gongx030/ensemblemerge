#' Get a neighbor graph
#' 
#' @param data a Seurat object
#' @param params a BaseMerge object
#' @param latent Whether or not use the latent representation to build the SNN (default: FALSE)
#' @param k.param	Defines k for the k-nearest neighbor algorithm (default: 20L)
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
	function(
		data, 
		params, 
		latent = FALSE,
		k.param = 20L,
		...
	) {

		if (latent){
			stopifnot(!is.null(x@reductions[[params@name]]))
			reduction <- params@name
			dims <- 1:ncol(x@reductions[[params@name]])
		}else{
			stopifnot(!is.null(x@reductions[[params@umap_name]]))
			reduction <- params@umap_name
			dims <- 1:ncol(x@reductions[[params@umap_name]])
		}

   	data <- FindNeighbors(
			data, 
			compute.SNN = TRUE, 
			reduction = reduction,
			dims = dims,
			k.param = k.param,
			graph.name = c(params@knn_name, params@snn_name),
			verbose = FALSE
		)

		data
	}
)

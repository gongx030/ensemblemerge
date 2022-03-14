#' Merge Seurat objects by the EnsembleMerge method
#'
#' @param x a Seurat object 
#' @param params a EnsembleMerge object
#' @param ... Additional arguments 
#' @return returns a Seurat object of the integrated data
#' @importFrom uwot umap
#' @importFrom Seurat DefaultAssay RunUMAP FindNeighbors
#' @importFrom methods as
#' @export
#'
setMethod(
	'ensemble',
	signature(
		x = 'Seurat',
		params = 'EnsembleMerge'
	),
	function(
		x,
		params,
		...
	){

		valid(x, params)

		n_reductions <- length(params@constituent_reduction_names)

		sprintf('ensemble | FindNeighbors | k.param=%d | latent=%s', params@k.param, params@latent) %>% message()
		nn <- list()
		for (i in 1:length(params@constituent_reduction_names)){

			m <- params@constituent_reduction_names[i]
			dims <- 1:ncol(x@reductions[[m]])

			x <- FindNeighbors(
				x,
				compute.SNN = TRUE,
				reduction = m,
				dims = dims,
				k.param = params@k.param,
				graph.name = c(params@constituent_knn_names[i], params@constituent_snn_names[i]),
				verbose = FALSE
			)

			nn[[i]] <- x[[params@constituent_snn_names[i]]]
		}

		if (n_reductions > 1){
			agreement <- Reduce('+', lapply(nn, function(x) x > 0)) 
			agreement <- as(agreement, 'dgCMatrix')
			agreement <- agreement / n_reductions
			wt <- sapply(1:n_reductions, function(i) sum(nn[[i]] * agreement) / sqrt(sum(nn[[i]] * nn[[i]]) * sum(agreement * agreement)))	
			wt <- (wt -  min(wt)) / (max(wt) - min(wt))
#			sigmoid <- function(x, a = 1, b = 1){1/(1+a*(x/(1-x))^-(b*2))}
#			wt <- sigmoid(wt)
			wt <- (wt -  min(wt)) / (max(wt) - min(wt))
			x@graphs[[params@snn_name]]<- Reduce('+', lapply(1:n_reductions, function(i) nn[[i]] * wt[i]))
			x@graphs[[params@snn_name]] <- x@graphs[[params@snn_name]] / sum(wt)
			names(wt) <- params@constituent_knn_names
			x@misc[['weight']][[params@snn_name]] <- wt
		}else{
			x@graphs[[params@snn_name]] <- as(nn[[1L]], 'dgCMatrix')
		}
	
		sprintf('ensemble | running UMAP') %>% message()
		x[[params@umap_name]] <- RunUMAP(
			x@graphs[[params@snn_name]],	# must be symmetric
			assay = params@raw_assay,
			n.components = params@umap_dim,
			seed.use = 1,
			verbose = FALSE
		)
		x
	}
)

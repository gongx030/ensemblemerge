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
#' @param data a SeuratList object 
#' @param params a SeuratMerge object 
#' @param ... Additional arguments
#' @return returns a data object with integrated data
#' @importFrom methods is
#'
setMethod(
	'Merge',
	signature(
		data  = 'SeuratList',
		params = 'SeuratMerge'
	),
	function(
		data,
		params,
		...
	){
   	data <- run_Seurat(params, data)
		data
	}
)

#' Merge Seurat objects
#'
#' @param data a SeuratList object 
#' @param params a SeuratMerge object 
#' @param ... Additional arguments
#' @return returns a data object with integrated data
#' @importFrom methods is
#'
setMethod(
	'Merge',
	signature(
		data  = 'Seurat',
		params = 'SeuratMerge'
	),
	function(
		data,
		params,
		...
	){
		stopifnot(valid(data, params))
		data <- SplitObject(data, split.by = params@batch)
   	data <- run_Seurat(params, data)
		data
	}
)

#' Merge Seurat objects
#'
#' @param data Seurat object containing single cell counts matrix
#' @param params a HarmonyMerge object 
#' @param ... Additional arguments
#' @return returns a Seurat object with integrated data
#'
setMethod(
	'Merge',
	signature(
		data  = 'Seurat',
		params = 'HarmonyMerge'
	),
	function(
		data,
		params,
		...
	){
		stopifnot(valid(data, params))
   	data <- run_Harmony(params, data)
		data
	}
)

#' Merge Seurat objects
#'
#' @param data Seurat object containing single cell counts matrix
#' @param params a FastMNNMerge object 
#' @param ... Additional arguments
#' @return returns a Seurat object with integrated data
#'
setMethod(
	'Merge',
	signature(
		data  = 'Seurat',
		params = 'FastMNNMerge'
	),
	function(
		data,
		params,
		...
	){
		stopifnot(valid(data, params))
   	data <- run_fastMNN(params, data)
		data
	}
)

#' Merge Seurat objects
#'
#' @param data Seurat object containing single cell counts matrix
#' @param params a UncorrectedMerge object 
#' @param ... Additional arguments
#' @return returns a Seurat object with integrated data
#'
setMethod(
	'Merge',
	signature(
		data  = 'Seurat',
		params = 'UncorrectedMerge'
	),
	function(
		data,
		params,
		...
	){
		stopifnot(valid(data, params))
   	data <- run_Uncorrected(params, data)
		data
	}
)

#' Merge Seurat objects
#'
#' @param data Seurat object containing single cell counts matrix
#' @param params a LigerMerge object 
#' @param ... Additional arguments
#' @return returns a Seurat object with integrated data
#'
setMethod(
	'Merge',
	signature(
		data  = 'Seurat',
		params = 'LigerMerge'
	),
	function(
		data,
		params,
		...
	){
		stopifnot(valid(data, params))
   	data <- run_Liger(params, data)
		data
	}
)

#' Merge Seurat objects
#'
#' @param data Seurat object containing single cell counts matrix
#' @param params a BBKNNMerge object 
#' @param ... Additional arguments
#' @return returns a Seurat object with integrated data
#'
setMethod(
	'Merge',
	signature(
		data  = 'Seurat',
		params = 'BBKNNMerge'
	),
	function(
		data,
		params,
		...
	){
		stopifnot(valid(data, params))
   	data <- run_BBKNN(params, data)
		data
	}
)

#' Merge Seurat objects
#'
#' @param data Seurat object containing single cell counts matrix
#' @param params a ScanoramaMerge object 
#' @param ... Additional arguments
#' @return returns a Seurat object with integrated data
#'
setMethod(
	'Merge',
	signature(
		data  = 'Seurat',
		params = 'ScanoramaMerge'
	),
	function(
		data,
		params,
		...
	){
		stopifnot(valid(data, params))
   	data <- run_Scanorama(params, data)
		data
	}
)

#' Merge Seurat objects
#'
#' @param data Seurat object containing single cell counts matrix
#' @param params a scVIMerge object 
#' @param ... Additional arguments
#' @return returns a Seurat object with integrated data
#'
setMethod(
	'Merge',
	signature(
		data  = 'Seurat',
		params = 'scVIMerge'
	),
	function(
		data,
		params,
		...
	){
		stopifnot(valid(data, params))
   	data <- run_scVI(params, data)
		data
	}
)

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
			d <- as.data.frame(rowData(data))
			seurat[[params@preprocess@raw_assay]]@meta.features <- d
		}

		Merge(seurat, params, ...)
})

#' Merge Seurat objects
#'
#' @param data a Seurat object 
#' @param params a MethodList object 
#' @param ... Additional arguments
#' @importFrom methods is
#' @importFrom uwot umap
#' @importFrom Seurat RunUMAP
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


		for (i in 1:length(params@constituent)){

			sprintf('Merge | Preprocess') %>% message()
			d <- Preprocess(data, params@preprocess, params@constituent[[i]])

			sprintf('Merge | running %s', params@constituent[[i]]@name) %>% message()
			d <- Merge(d, params@constituent[[i]])

			# it is likely that ncol(d) will be smaller than ncol(data) due to the cell-wise filtering
			# in this case, only the cells present in `d` will be kept in `data`
			if (ncol(d) < ncol(data))
				data <- data[, colnames(d)]

			if (is(params@constituent[[i]], 'BBKNNMerge')){

				# BBKNN does not have latent representations
				data[[params@constituent[[i]]@umap_name]] <- d[[params@constituent[[i]]@umap_name]]

			}else{ # where the latent representations are available

				data[[params@constituent[[i]]@name]] <- CreateDimReducObject(
					embeddings = d@reductions[[params@constituent[[i]]@name]]@cell.embeddings[colnames(data), , drop = FALSE],
					assay = params@constituent[[i]]@raw_assay,
					key = params@constituent[[i]]@reduction_key
				)

				sprintf('Merge | running UMAP') %>% message()
				data <- RunUMAP(
					data, 
					reduction = params@constituent[[i]]@name, 
					dims = 1:params@constituent[[i]]@npcs, 
					n.components = params@constituent[[i]]@umap_dim,
					reduction.key = params@constituent[[i]]@umap_key, 
					reduction.name = params@constituent[[i]]@umap_name, 
					seed.use = 1, 
					verbose = FALSE
				)
			}
		}
		data
	}
)
	

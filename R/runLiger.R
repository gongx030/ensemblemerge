#' Run Liger merging
#'
#' @param params a LigerParams object
#' @param data a Seurat object
#'
#' @return returns a Seurat object of the integrated data
#'
#' @importFrom Seurat GetAssayData DefaultAssay SplitObject CreateDimReducObject
#'
run_Liger <- function(params, data){

  assay <- DefaultAssay(data)

	d <- SplitObject(data, split.by = params@batch)
	counts <- lapply(d, function(x){
		GetAssayData(
      object = x,
			slot = 'counts',
			assay = assay
		)
	})

	# Here we used the rliger pipeline for noramlization and scaling
	# Using scaled data from the univeral pipeline (e.g. Seurat) will sometimes have error 
	# see https://github.com/welch-lab/liger/issues/174
	# However, we still use the HVGs selected by the univeral pipeline
	#
	d <- rliger::createLiger(counts, remove.missing = FALSE)
	d <- rliger::normalize(d, remove.missing = FALSE)
	d@var.genes <- VariableFeatures(data)
	d <- rliger::scaleNotCenter(d)
	d <- rliger::optimizeALS(d, k = params@npcs, verbose = FALSE)
	d <- rliger::quantile_norm(d, verbose = FALSE)

	data[[params@name]] <- CreateDimReducObject(
		embeddings = d@H.norm,
		assay = assay,
		key = params@reduction_key
	)
	data
}

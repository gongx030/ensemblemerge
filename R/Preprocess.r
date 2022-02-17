#' Preprocess a Seurat objects by the Seurat pipeline
#'
#' Adopted from https://satijalab.org/seurat/articles/integration_introduction.html
#'
#' @param params a SeuratPreprocess object 
#' @param data a Seurat object
#' @param ... Additional arguments
#' @return returns a Seurat object

#' @importFrom Seurat SplitObject NormalizeData FindVariableFeatures ScaleData SelectIntegrationFeatures
#'
setMethod(
	'Preprocess',
	signature(
		params = 'SeuratPreprocess',
		data  = 'Seurat'
	),
	function(
		params,
		data,
		...
	){

		cn <- colnames(data)

		stopifnot(!any(duplicated(cn)))

		if (params@batchwise){
			batch_list <- SplitObject(data, split.by = params@batch)
		}else{
			batch_list <- list(data)
		}

	  for(i in 1:length(batch_list)) {
			if(params@norm_data){
				batch_list[[i]] <- NormalizeData(
					batch_list[[i]],
					normalization.method = params@norm_method,
					scale.factor = params@scale_factor,
					verbose = FALSE
				)
	    }
	    batch_list[[i]] <- FindVariableFeatures(
				batch_list[[i]],
				selection.method = params@selection.method,
				nfeatures = params@numHVG,
				verbose = FALSE
			)
		}

		# select features that are repeatedly variable across datasets for integration
		features <- SelectIntegrationFeatures(batch_list, verbose = FALSE)

		data@assays[[data@active.assay]]@var.features <- features

		if(params@regressUMI && params@scaling) {
			data <- ScaleData(object = data, vars.to.regress = params@vars.to.regress, verbose = FALSE)
		} else if (params@scaling) { # default option
			data <- ScaleData(object = data, verbose = FALSE)
		}
		data

	}
)


#' Preprocess a Seurat object by the Scanpy pipeline
#'
#' @param params a ScanpyPreprocess object
#' @param data  a Seurat object
#' @param ... Additional arguments
#' @return returns a Seurat object
#'
setMethod(
	'Preprocess',
	signature(
		params = 'ScanpyPreprocess',
		data  = 'Seurat'
	),
	function(
		params,
		data,
		...
	){

		data <- as.SingleCellExperiment(data)
		Preprocess(params, data, ...)
	}
)

#' Preprocess a SingleCellExperiment object by the Scanpy pipeline
#'
#' @param params a ScanpyPreprocess object
#' @param data a SingleCellExperiment object
#' @return returns a SingleCellExperiment object
#' @param ... Additional arguments
#' @importFrom reticulate import
#'
setMethod(
	'Preprocess',
	signature(
		params = 'ScanpyPreprocess',
		data  = 'SingleCellExperiment'
	),
	function(
		params,
		data,
		...
	){

		h5ad_file <- tempfile(fileext = '.h5ad')
		data = sceasy::convertFormat(data, from = "sce", to = "anndata", out = h5ad_file)
		sc <- import("scanpy")
		sc$pp$filter_cells(data, min_genes = params@min_genes)
		sc$pp$filter_genes(data, min_cells = params@min_cells)
		sc$pp$normalize_total(data, target_sum = params@scale_factor)
		sc$pp$log1p(data)
		sc$pp$highly_variable_genes(data, min_mean = params@min_mean, max_mean= params@max_mean, min_disp = params@min_disp)
		sc$pp$scale(data)
		sc$tl$pca(data)
		data$write(filename = h5ad_file)
		data <- zellkonverter::readH5AD(h5ad_file, reader = 'R')
		rowData(data)$highly_variable <- as.logical(rowData(data)$highly_variable)
		unlink(h5ad_file)
		data
	}
)

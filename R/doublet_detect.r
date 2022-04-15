setClass(
	'scDblFinderDoubletDetect',
	representation(
		dbr = 'numeric',
		nfeatures = 'integer',
		dims = 'integer',
		includePCs = 'integer',
		score = 'character',
		metric = 'character',
		nrounds = 'numeric',
		max_depth = 'integer',
		iter = 'integer'
	),
	contains = c('BaseDoubletDetect'),
	prototype(
		dbr = 0.1,
		nfeatures = 1000L,
		dims = 20L,
		includePCs = 10L,
		score = 'xgb',
		metric = 'logloss',
		nrounds = 0.25,
		max_depth = 4L,
		iter = 3L,
		dependences = list(
			new('RPackage', package_name = 'scDblFinder', package_version = '1.8.0')
		)
	),
)

#' Doublet detection by scDblFinder (https://bioconductor.org/packages/release/bioc/html/scDblFinder.html)
#'
#' @param x a Seurat object
#' @param params a scDblFinderDoubletDetectobject
#' @param ... Additional arguments
#' @return returns a data object with doublet removed
#' @references Germain P, Lun A, Macnair W, Robinson M (2021). “Doublet identification in single-cell sequencing data using scDblFinder.” f1000research. doi: 10.12688/f1000research.73600.1.
#'
setMethod(
	'DetectDoublet',
	signature(
		x = 'Seurat',
		params = 'scDblFinderDoubletDetect'
	),
	function(
		x,
		params,
		...
	){

		raw_assay <- params@raw_assay
		sce <- SingleCellExperiment(
			assays = list(
				counts = x@assays[[raw_assay]]@counts 
			)
		)

		sce <- scDblFinder::scDblFinder(
			sce,
			clusters = NULL,
			samples = NULL,
			clustCor = NULL,
			artificialDoublets = NULL,
			knownDoublets = NULL,
			knownUse = 'discard',
			dbr = params@dbr,
			dbr.sd = NULL,
			nfeatures = params@nfeatures,
			dims = params@dims,
			k = NULL,
			removeUnidentifiable = TRUE,
			includePCs = params@includePCs,
			propRandom = 0,
			propMarkers = 0,
			aggregateFeatures = FALSE,
			returnType = 'sce',
			score = params@score,
			processing = "default",
			metric = params@metric,
			nrounds = params@nrounds,
			max_depth = params@max_depth,
			iter = params@iter,
			trainingFeatures = NULL,
			multiSampleMode = 'asOne',
			threshold = TRUE,
			verbose = FALSE
		)

		is_singlet <- colData(sce)$scDblFinder.class == 'singlet'
		sprintf('DetectDoublet | removing %d doublets', sum(!is_singlet)) %>% message()
		x <- x[, is_singlet]
		x
	}
)


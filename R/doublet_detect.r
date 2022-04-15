#' Remove doublets from a SeuratList object
#'
#' @param x a SeuratList object
#' @param params a BaseDoubletDetect object
#' @param ... Additional arguments
#' @return returns a SeuratList object with doublets removed
#'
setMethod(
	'DetectDoublet',
	signature(
		x = 'SeuratList',
		params = 'BaseDoubletDetect'
	),
	function(
		x,
		params,
		...
	){

		for (i in 1:length(x)){
			x[[i]] <- DetectDoublet(x[[i]], params, ...)
		}
		x
	}
)


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

		raw_assay <- params@preprocess@raw_assay
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


setClass(
	'DoubletFinderoubletDetect',
	representation(
		dbr = 'numeric',
		ndims = 'integer',
		sct = 'logical',
		pN = 'numeric'
	),
	contains = c('BaseDoubletDetect'),
	prototype(
		dbr = 0.1,
		ndims = 10L,
		pN = 0.25,
		dependences = list(
			new('RPackage', package_name = 'DoubletFinder', package_version = '2.0.3')
		)
	),
)

#' @importFrom methods callNextMethod is 
#'
setMethod('initialize', 'DoubletFinderoubletDetect', function(.Object, check_dependencies = TRUE, ...){
	.Object <- callNextMethod(.Object, check_dependencies = check_dependencies, ...)
	if (is(.Object@preprocess, 'SeuratPreprocess')){
		.Object@sct <- FALSE
	}else{
		.NotYetImplemented()
	}
	.Object
})

#' Doublet detection by DoubletFinderoubletDetect (https://github.com/chris-mcginnis-ucsf/DoubletFinder)
#'
#' @param x a Seurat object
#' @param params a DoubletFinderoubletDetect 
#' @param ... Additional arguments
#' @return returns a data object with doublet removed
#' @references McGinnis CS, Murrow LM, Gartner ZJ. DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors. Cell Syst. 2019 Apr 24;8(4):329-337.e4. doi: 10.1016/j.cels.2019.03.003. Epub 2019 Apr 3. PMID: 30954475; PMCID: PMC6853612.
#' @importFrom grDevices png dev.off
#'

setMethod(
	'DetectDoublet',
	signature(
		x = 'Seurat',
		params = 'DoubletFinderoubletDetect'
	),
	function(
		x,
		params,
		...
	){

		params_pca <- new('PCAEmbed', preprocess = params@preprocess, ndims = params@ndims)
		x2 <- Embed(x, params_pca)
		res <- DoubletFinder::paramSweep_v3(x2, PCs = 1:params@ndims, sct = params@sct)
		res_sweep <- DoubletFinder::summarizeSweep(res, GT = FALSE)

		# suppress the plotting (https://stackoverflow.com/questions/20363266/how-can-i-suppress-the-creation-of-a-plot-while-calling-a-function-in-r)
		ff <- tempfile()
		png(filename=ff)
		res_pk <- DoubletFinder::find.pK(res_sweep)
		dev.off()
		unlink(ff)

		pK <- res_pk$pK[which.max(res_pk$BCmetric)] %>% as.character() %>% as.numeric()

		x2 <- DoubletFinder::doubletFinder_v3(
			x2,
			PCs = 1:params@ndims, 
			pN = params@pN, 
			pK = pK,
			nExp = round(params@dbr * ncol(x2)),
			reuse.pANN = FALSE, 
			sct = params@sct
		)

		is_singlet <- x2@meta.data[[grep('^DF.classifications', colnames(x2@meta.data))[1]]] == 'Singlet'
		sprintf('DetectDoublet | removing %d doublets', sum(!is_singlet)) %>% message()
		x <- x[, is_singlet]
		x

	}
)


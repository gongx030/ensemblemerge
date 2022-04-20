#' Remove ambient RNAs from a SeuratList object
#'
#' @param x a SeuratList object
#' @param params a BaseAmbientRNARemoval object
#' @param ... Additional arguments
#' @return returns a SeuratList object with ambient RNA removed
#' @export
#'
setMethod(
	'RemoveAmbientRNA',
	signature(
		x = 'SeuratList',
		params = 'BaseAmbientRNARemoval'
	),
	function(
		x,
		params,
		...
	){

		for (i in 1:length(x)){
			x[[i]] <- RemoveAmbientRNA(x[[i]], params, ...)
		}
		x
	}
)

setClass(
	'decontXRemoveAmbientRNA',
	representation(
	),
	contains = c('BaseAmbientRNARemoval'),
	prototype(
		name = 'decontXRemoveAmbientRNA',
		dependences = list(
			new('RPackage', package_name = 'celda', package_version = '1.10.0')
		)
	)
)

#' Remove ambient RNAs by decontX (http://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html)
#'
#' @param x a Seurat object
#' @param params a decontXRemoveAmbientRNA object
#' @param ... Additional arguments
#' @return returns a data object with ambient RNA removed
#' @references Yang, S., Corbett, S.E., Koga, Y. et al. Decontamination of ambient RNA in single-cell RNA-seq with DecontX. Genome Biol 21, 57 (2020). https://doi.org/10.1186/s13059-020-1950-6
#' @export
#'
setMethod(
	'RemoveAmbientRNA',
	signature(
		x = 'Seurat',
		params = 'decontXRemoveAmbientRNA'
	),
	function(
		x,
		params,
		...
	){

		raw_assay <- params@normalize@assay_name

		sce <- SingleCellExperiment(
			assays = list(
				counts = x@assays[[raw_assay]]@counts 
			)
		)

		params_embed <- new('PCAEmbed', normalize = params@normalize)
		x <- Embed(x, params_embed)

	  params_cluster <- new('LouvainCluster', embedding = params_embed)
		x <- Cluster(x, params_cluster)

		sce <- celda::decontX(
			sce,
			assayName = 'counts',
			z = x[[params_cluster@cluster_name]][, 1],
			batch = NULL,
			background = NULL,
			bgAssayName = NULL,
			maxIter = 500,
			delta = c(10, 10),
			estimateDelta = TRUE,
			convergence = 0.001,
			iterLogLik = 10,
			seed = params@seed,
			logfile = NULL,
			verbose = FALSE 
		)

		x[[params@normalize@assay_name]] <- CreateAssayObject(
			counts = assays(sce)$decontXcounts
		)

		x <- Normalize(x, params@normalize)
		x
	}
)


#' Annotate a SeuratList
#'
#' @param x a SeuratList object
#' @param params a BaseAnnotate object
#' @param ... Additional arguments
#' @return returns a SeuratList object with cell-wise annotation
#' @export
#'
setMethod(
	'Annotate',
	signature(
		x = 'SeuratList',
		params = 'BaseAnnotate'
	),
	function(
		x,
		params,
		...
	){

		for (i in 1:length(x)){
			x[[i]] <- Annotate(x[[i]], params, ...)
		}
		x
	}
)


#' The scCATCHAnnotate class
#'
setClass(
	'scCATCHAnnotate',
	representation(
		gene_marker = 'BaseGeneMarkers',
		cluster = 'BaseCluster',
		tissue = 'character',
		species = 'character',
		cancer = 'character',
		cell_min_pct = 'numeric',
		logfc = 'numeric',
		pvalue = 'numeric'

	),
	contains = c('BaseAnnotate'),
	prototype(
		name = 'scCATCHAnnotate',
		cancer = 'Normal',
		dependences = list(
			new('RPackage', package_name = 'scCATCH', package_version = '3.0')
		),
		cell_min_pct = 0.25,
		logfc = 0.25,
		pvalue = 0.05
	)
)


#' Annotate a Seurat object with scCATCH (https://cran.r-project.org/web/packages/scCATCH/index.html)
#'
#' @param x a Seurat object
#' @param params a scCATCHAnnotate object
#' @param ... Additional arguments
#' @return returns a data object with cell-wise annotation
#' @importFrom rlang .data
#' @export
#'
setMethod(
	'Annotate',
	signature(
		x = 'Seurat',
		params = 'scCATCHAnnotate'
	),
	function(
		x,
		params,
		...
	){

		# yet to be implemented
		# stopifnot(valid(x, params))

		cellmatch  <- data.frame(
			species = NA,
			tissue = NA,
			cancer = NA,
			conditoin = NA,
			subtype1 = NA,
			subtype2 = NA,
			subtype3 = NA,
			celltype = params@gene_marker@celltype[[params@gene_marker@level]],
			gene = params@gene_marker@celltype$gene,
			resource = params@gene_marker@name,
			pmid = NA
		)
		raw_assay <- params@cluster@embedding@normalize@assay_name
		cls <- x[[params@cluster@cluster_name]][, 1] %>%
			as.character()

		obj <- scCATCH::createscCATCH(data = x@assays[[raw_assay]]@data, cluster = cls)
		obj <- scCATCH::findmarkergene(
			object = obj, 
			if_use_custom_marker = TRUE,
			marker = cellmatch,
			cell_min_pct = params@cell_min_pct,
			logfc = params@logfc,
			pvalue = params@pvalue,
			verbose = FALSE
		)
		obj <- scCATCH::findcelltype(object = obj, verbose = FALSE)

		rownames(obj@celltype) <- obj@celltype[, 'cluster']
		x@meta.data[[params@annotate_name]] <- obj@celltype[obj@meta$cluster, 'cell_type']
		x
	}
)

#' The SCINAAnnotate class
#'
setClass(
	'SCINAAnnotate',
	representation(
		species = 'character',
		cancer = 'character',
		max_iter = 'integer',
		convergence_n = 'integer',
		convergence_rate = 'numeric',
		sensitivity_cutoff = 'numeric',
		rm_overlap = 'logical',
		allow_unknown = 'logical'
	),
	contains = c('BaseAnnotate'),
	prototype(
		name = 'SCINAAnnotate',
		cancer = 'Normal',
		max_iter = 100L,
		convergence_n = 10L,
		convergence_rate = 0.999,
		sensitivity_cutoff = 0.9,
		rm_overlap = FALSE,
		allow_unknown = TRUE,
		dependences = list(
			new('RPackage', package_name = 'SCINA', package_version = '1.2.0'),
			new('RPackage', package_name = 'scCATCH', package_version = '3.0')
		)
	),
)

setMethod('initialize', 'SCINAAnnotate', function(.Object, check_dependencies = TRUE, ...){
	.Object <- callNextMethod(.Object, check_dependencies = check_dependencies, ...)
	if (.Object@genome == 'mm10'){
		.Object@species <- 'Mouse'
	}else if (.Object@genome == 'hg19'){
		.Object@species <- 'Human'
	}
	.Object
})


#' Annotate a Seurat object with SCINA (https://github.com/jcao89757/SCINA)
#'
#' @param x a Seurat object
#' @param params a SCINAAnnotate object
#' @param ... Additional arguments
#' @return returns a data object with cell-wise annotation
#' @importFrom rlang .data
#' @export
#'
setMethod(
	'Annotate',
	signature(
		x = 'Seurat',
		params = 'SCINAAnnotate'
	),
	function(
		x,
		params,
		...
	){

		# yet to be implemented
		# stopifnot(valid(x, params))

		signatures <- scCATCH::cellmatch %>%
			filter(.data$cancer == params@cancer & .data$species == params@species) %>%
			select(.data$gene, .data$celltype)

		signatures <- split(signatures$gene, list(signatures$celltype))

		raw_assay <- params@normalize@assay_name

		results <- SCINA::SCINA(
			x@assays[[raw_assay]]@data, 
			signatures,
			max_iter = params@max_iter,
			convergence_n = params@convergence_n,
			convergence_rate = params@convergence_rate,
			sensitivity_cutoff = params@sensitivity_cutoff,
			rm_overlap = params@rm_overlap,
			allow_unknown = params@allow_unknown
		)

		x@meta.data[[params@annotate_name]] <- results$cell_labels
		x
	}
)

setClass(
	'CellAssignAnnotate',
	representation(
		gene_marker = 'BaseGeneMarkers'
	),
	contains = c('BaseAnnotate'),
	prototype(
		name = 'CellAssignAnnotate',
		dependences = list(
			new('RPackage', package_name = 'tensorflow', package_version = '2.8.0'),
			new('PythonPackage', package_name = 'tensorflow', package_version = '2.8.0'),
			new('RPackage', package_name = 'cellassign', package_version = '0.99.21'),
			new('RPackage', package_name = 'scran', package_version = '1.22.1'),
			new('PythonPackage', package_name = 'tensorflow-probability', package_version = '0.15.0'),
			new('RPackage', package_name = 'tfprobability', package_version = '0.15.0')
		)
	),
)

#' Annotate a Seurat object with CellAssignAnnotate (https://github.com/Irrationone/cellassign)
#'
#' @param x a Seurat object
#' @param params a CellAssignAnnotate object
#' @param ... Additional arguments
#' @return returns a data object with cell-wise annotation
#' @references Zhang, A.W., O’Flanagan, C., Chavez, E.A. et al. Probabilistic cell-type assignment of single-cell RNA-seq for tumor microenvironment profiling. Nat Methods 16, 1007–1015 (2019). https://doi.org/10.1038/s41592-019-0529-1
#' @importFrom Seurat GetAssayData
#'
setMethod(
	'Annotate',
	signature(
		x = 'Seurat',
		params = 'CellAssignAnnotate'
	),
	function(
		x,
		params,
		...
	){

		.NotYetImplemented()

		# it is still possible that there are genes not detected in any cells in current batch since the initial filtering
		# was done on the full datasets, before splitting based upon batch indicator
#		rc <- rowSums(GetAssayData(x, 'data') > 0)
#		invalid <- rc < params@normalize@preprocess@min_cells
#		if (any(invalid)){
#			x <- x[!invalid, ]
#		}

#		raw_assay <- params@normalize@assay_name

		# it appears that SC3 only takes dense matrix as the input
#		sce <- SingleCellExperiment(
#			assays = list(
#				counts = x@assays[[raw_assay]]@counts 
#			)
#		)
#		sce <- scran::computeSumFactors(sce)
#		s <- sizeFactors(sce)

#		gs <- intersect(rownames(x), params@gene_marker@celltype[, 'gene'])
#		gm <- params@gene_marker@celltype[params@gene_marker@celltype[, 'gene'] %in% gs, ]
#		cs <- gm[, params@gene_marker@level]
#		mm <- sparseMatrix(
#			i = factor(gm[, 'gene'], gs) %>% as.numeric(),
#			j = factor(cs, unique(cs)) %>% as.numeric(),
#			dims = c(length(gs), length(unique(cs))),
#			dimnames = list(gs, unique(cs))
#		) %>% 
#			as.matrix()

#		sce <- sce[gs, ]

#		results <- cellassign::cellassign(
#			exprs_obj = sce,
#			marker_gene_info = mm,
#			s = s, 
#			learning_rate = 1e-2, 
#			shrinkage = TRUE,
#			verbose = FALSE
#		)

	}
)


setClass(
	'scMRMAAnnotate',
	representation(
		cluster = 'BaseCluster',
		gene_marker = 'PanglaoDBGeneMarkers'
	),
	contains = c('BaseAnnotate'),
	prototype(
		name = 'scMRMAAnnotate',
		dependences = list(
			new('RPackage', package_name = 'scMRMA', package_version = '1.0')
		)
	),
)

#' Annotate a Seurat object with scMRMAAnnotate (https://github.com/JiaLiVUMC/scMRMA)
#'
#' @param x a Seurat object
#' @param params a CellAssignAnnotate object
#' @param ... Additional arguments
#' @return returns a data object with cell-wise annotation
#' @references Jia Li, Quanhu Sheng, Yu Shyr, Qi Liu, scMRMA: single cell multiresolution marker-based annotation, Nucleic Acids Research, Volume 50, Issue 2, 25 January 2022, Page e7, https://doi.org/10.1093/nar/gkab931
#'
setMethod(
	'Annotate',
	signature(
		x = 'Seurat',
		params = 'scMRMAAnnotate'
	),
	function(
		x,
		params,
		...
	){

		.NotYetImplemented()

		# need to make sure that the rownames of matrix are gene symbols

#		raw_assay <- params@normalize@assay_name
#		h <- x@assays[[raw_assay]]@meta.features[[params@normalize@feature_field]]
#	 	cls <- x[[params@cluster@cluster_name]][, 1]
#		names(cls) <- colnames(x)

#		scMRMA::scMRMA(
#			x@assays[[raw_assay]]@data[h, ],
#			species = NULL,
#			db = NULL,
#			selfDB = params@gene_marker@celltype,
#			p = 0.05,
#			normalizedData = TRUE,
#			selfClusters = cls,
#			k = 20
#		)

		# following error when running `scMRMA`
		# User-provided cell type database will be used.
		# Multi Resolution Annotation Started.
		# Level 1 annotation started.
		# Error in matrix(NA, nrow = nrow(sub), ncol = length(unique(clusters))) :
		# non-numeric matrix extent

	}
)


setClass(
	'clustifyrAnnotate',
	representation(
		cluster = 'BaseCluster',
		gene_marker = 'PanglaoDBGeneMarkers',
		metric = 'character'
	),
	contains = c('BaseAnnotate'),
	prototype(
		name = 'clustifyrAnnotate',
		metric = 'jaccard',
		dependences = list(
			new('RPackage', package_name = 'clustifyr', package_version = '1.6.0')
		)
	)
)

#' Annotate a Seurat object with clustifyrAnnotate (https://github.com/rnabioco/clustifyr)
#'
#' @param x a Seurat object
#' @param params a clustifyrAnnotate object
#' @param ... Additional arguments
#' @return returns a data object with cell-wise annotation
#' @references Fu R, Gillen AE, Sheridan RM, et al. clustifyr: an R package for automated single-cell RNA sequencing cluster classification. F1000Res. 2020;9:223. Published 2020 Apr 1. doi:10.12688/f1000research.22969.2
#' @export
#'
setMethod(
	'Annotate',
	signature(
		x = 'Seurat',
		params = 'clustifyrAnnotate'
	),
	function(
		x,
		params,
		...
	){

		raw_assay <- params@normalize@assay_name
		features <- x@assays[[raw_assay]]@var.features

		marker <- split(params@gene_marker@celltype[, 'gene'], list(params@gene_marker@celltype[, params@gene_marker@level]))

		cor_mat <- clustifyr::clustify_lists(
			input = x@assays[[raw_assay]]@data[features, ],	# requires normalized counts
			metadata = x@meta.data,            # meta.data table containing cell clusters
			cluster_col = params@cluster@cluster_name,
			marker = marker,                 # list of known marker genes
			metric = params@metric,                   # test to use for assigning cell types
			verbose = FALSE
		)

		results <- clustifyr::cor_to_call(
			cor_mat,
			metadata =  x@meta.data,
			cluster_col = params@cluster@cluster_name,
			collapse_to_cluster = FALSE,
			threshold = 0
		)

		celltype <- results[['type']] %>% as.vector()
		names(celltype) <- results[[params@cluster@cluster_name]]
		x@meta.data[[params@annotate_name]] <- celltype[as.character(x@meta.data[[params@cluster@cluster_name]])]
		x

	}
)

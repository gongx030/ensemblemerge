#' parent
#' 
#' set base params object
#'
#' @slot batch character name of batch in dataset metadata
#' @slot dimreduc_names name of dimension reduction in metadata, i.e. "PCA"
#' @slot return data type to return, can be Seurat or SingleCellExperiment 
#'
#' @export
setClass(
	'parent', 
	representation(
		batch = "character",
    dimreduc_names = "character",
    return = "character",
    name = "character"
	),
  prototype(batch = "batch",
            dimreduc_names = c("PCA" = "pca",
                              "UMAP" = "umap",
                              "tSNE" = "tsne"),
            return = "SingleCellExperiment",
            name = "parent")
)

#' SeuratMerge
#'
#' class for merging dataset with Seurat
#'
#' @slot npcs number of principle components to use in dimension reduction
#' @slot seed value to set for seed for umap
#' @slot dims number of dimensions used in seurat preprocessing
#' @slot k.weight weight for neighbor function
#'
#' @export
setClass(
	'SeuratMerge', 
	representation(
		npcs = "numeric",
    seed = "numeric",
    dims = "numeric",
    k.weight = "numeric"
	),
	contains = c('parent'),
  prototype(npcs = 20, 
            seed = 10,
            dims = 20,
            k.weight = 100)
)

#' SeuratNormalize
#'
#' normalization class for Seurat
#'
#' @slot norm_data boolean check to run data normalization
#' @slot scaling boolean check to scale data
#' @slot regressUMI boolean check to regress the UMI counts
#' @slot min_cells threshold to set the minimum number of cells a gene needs to be present in
#' @slot min_genes threshold of minimum number of genes in cells
#' @slot norm_method set normalization method, default is "LogNormalize"
#' @slot scale_factor factor for scaling
#'
#' @export
setClass(
	'SeuratNormalize', 
	representation(
		norm_data = "logical",
    scaling = "logical",
    regressUMI = "logical",
    min_cells = "numeric",
    min_genes = "numeric",
    norm_method = "character",
    scale_factor = "numeric"
	),
	contains = c('parent'),
  prototype(norm_data = TRUE,
            scaling = TRUE,
            regressUMI = FALSE,
            min_cells = 10,
            min_genes = 300,
            norm_method = "LogNormalize",
            scale_factor = 10000)
)

#' SeuratHVG
#'
#' Highly variable gene selection class for Seurat
#'
#' @slot numVG number of variable genes for integration
#' @slot numHVG number of highly variable genes to select for downstream analysis
#'
#'
#' @export
setClass(
	'SeuratHVG', 
	representation(
    numVG = "numeric",
    numHVG = "numeric"
	),
	contains = c('parent'),
  prototype(numVG = 300,
            numHVG = 2000)
)

#' SeuratParams
#'
#' base class for Seurat processing
#'
#' @slot altExp_names name(s) for additional experiment object in dataset to keep in downstream analysis
#'
#'
#' @export
setClass(
	'SeuratParams', 
	representation(altExp_names = "character"),
	contains = c('SeuratNormalize',
                'SeuratHVG',
                'SeuratMerge'),
  prototype(altExp_names = "RNA",
  name = "Seurat")
)

#' HarmonyMerge
#'
#' merge class for Harmony
#'
#' @slot theta_harmony diversity clustering penalty parameter, larger values increase diversity
#' @slot npcs number of principle components to use in umap
#' @slot seed set seed value for umap clustering
#' @slot num_clust number of clusters to use in harmony integration
#' @slot max_iter_cluster maximum number of learning iterations per cluster
#'
#'
#' @export
setClass(
	'HarmonyMerge', 
	representation(
    theta_harmony = "numeric",
    npcs = "numeric",
    seed = "numeric",
    num_clust = "numeric",
    max_iter_cluster = "numeric"
	),
	contains = c('parent'),
  prototype(theta_harmony = 2,
            num_clust = 50,
            max_iter_cluster = 100,
            npcs = 20,
            seed = 10)
)

#' HarmonyParams
#'
#' base parameters class for harmony integration
#'
#'
#' @export
setClass(
	'HarmonyParams', 
	representation(),
	contains = c('SeuratNormalize',
                'SeuratHVG',
                'HarmonyMerge'),
                prototype(name = "Harmony")
)

#' UncorrectedParams
#'
#' parameters class for uncorrected integration
#'
#' @slot vars_to_regress list of names or name of metadata elements to regress
#' @slot hvg boolean to check if hvg selection should be performed
#'
#'
#' @export
setClass(
	'UncorrectedParams', 
	representation(vars_to_regress = "character",
                hvg = "logical"),
	contains = c('SeuratNormalize',
                'SeuratHVG',
                'SeuratMerge'),
  prototype(vars_to_regress = c("nUMI"),
                                hvg = TRUE,
                                name = "Uncorrected")
)

#' FastMNNParams
#'
#' parameters class for fastMNN integration
#'
#'
#'
#' @export
setClass(
	'FastMNNParams', 
	representation(),
	contains = c('UncorrectedParams'),
  prototype(name = "FastMNN")
)

#' LigerNormalize
#'
#' parameters class for liger normalization
#'
#'
#' @export
setClass(
	'LigerNormalize', 
	representation(
	),
	contains = c('parent')
)

#' LigerHVG
#'
#' parameters class for liger highly variable gene
#'
#' @slot var_threshold theshold for variance of each gene in selecting variable features
#'
#' @export
setClass(
	'LigerHVG', 
	representation(
    var_threshold = "numeric"
	),
	contains = c('parent'),
  prototype(var_threshold = 0.1)
)

#' LigerMerge
#'
#' parameters class for liger merge
#'
#' @slot k number of clustrs for liger merge
#' @slot nrep number of restarts to perform
#' @slot lambda Regularization parameter. Larger values penalize dataset-specific effects more strongly
#'
#'
#' @export
setClass(
	'LigerMerge', 
	representation(
    k = "numeric",
    nrep = "numeric",
    lambda = "numeric"
	),
	contains = c('parent'),
  prototype(k = 20,
            nrep = 3,
            lambda = 5)
)

#' LigerParams
#'
#' parameters class for liger operations
#'
#' @export
setClass(
	'LigerParams', 
	representation(),
	contains = c('LigerNormalize',
                'LigerHVG',
                'LigerMerge',
                'UncorrectedParams'),
  prototype(name = "Liger")
)

#' BBKNNNormalize
#'
#' @export
setClass(
	'BBKNNNormalize', 
	representation(
    min_genes = "numeric",
    min_cells = "numeric",
    svd_solver = "character",
    scale_factor = "numeric",
    npcs = "numeric",
    nhvg = "numeric",
    n_neighbors = "numeric",
    confounder_key = "character"
	),
	contains = c('parent'),
  prototype(min_genes = 300,
            min_genes = 5,
            svd_solver = "arpack",
            scale_factor = 10000,
            npcs = 20,
            nhvg = 2000,
            n_neighbors = 10,
            confounder_key = "leiden")
)

#' BBKNNMerge
#'
#' @export
setClass(
	'BBKNNMerge', 
	representation(
    save_knn = "logical",
    copy = "logical",
    neighbors_within_batch = "numeric",
    approx = "logical",
    trim = "numeric",
    graph_name = "character",
    nn_name = "character",
    ridge_regress = "logical"
	),
	contains = c('parent'),
  prototype(save_knn = TRUE,
            copy = TRUE,
            neighbors_within_batch = 5,
            approx = FALSE,
            trim = 50,
            graph_name = "bbknn_graph",
            nn_name = "bbknn",
            ridge_regress = TRUE)
)

#' BBKNNParams
#'
#' @export
setClass(
	'BBKNNParams', 
	representation(),
	contains = c('BBKNNNormalize',
              'SeuratNormalize',
              'BBKNNMerge'),
  prototype(name = "bbknn")
)

#' ScanoramaParams
#'
#' @export
setClass(
	'ScanoramaParams', 
	representation(min_cells = "numeric",
          min_genes = "numeric",
          batch_size = "numeric",
          return_dense = "logical",
          knn = "numeric",
          svd_solver = "character",
          npcs = "numeric",
          nhvg = "numeric"),
	contains = c('parent'),
  prototype(min_genes = 300,
            min_cells = 5,
            batch_size = 30,
            return_dense = TRUE,
            knn = 10,
            svd_solver = "arpack",
            npcs = 20,
            nhvg = 2000,
            name = "Scanorama")
)


#' scVIParams
#'
#' @export
setClass(
	'scVIParams', 
	representation(min_cells = "numeric",
          min_genes = "numeric",
          batch_size = "numeric",
          return_dense = "logical",
          knn = "numeric",
          svd_solver = "character",
          npcs = "numeric",
          nhvg = "numeric",
          use_cuda = "logical"),
	contains = c('parent'),
  prototype(min_genes = 300,
            min_cells = 5,
            batch_size = 30,
            return_dense = TRUE,
            knn = 10,
            svd_solver = "arpack",
            npcs = 20,
            nhvg = 2000,
            name = "Scanorama",
            use_cuda = TRUE)
)
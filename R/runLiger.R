#' Run Liger merging
#'
#' @param params a LigerParams object
#' @param data a Seurat object
#'
#' @return returns a Seurat object of the integrated data
#'
run_Liger <- function(params, data){

  data <- data %>% 
		RunOptimizeALS(
			k = params@npcs, 
			lambda = params@lambda, 
			split.by = params@batch, 
			reduction.name = params@name, 
			reduction.key = params@name
		) %>%
 	 	RunQuantileNorm(
			split.by = params@batch, 
			reduction = params@name, 
			reduction.name = params@name, 
			reduction.key = params@name
		)
	data
}

#' Run optimizeALS on a Seurat object
#'
#' @param object A merged Seurat object
#' @param k see rliger::optimizeALS
#' @param assay  see rliger::optimizeALS
#' @param split.by Attribute for splitting, defaults to "orig.ident"
#' @param lambda see rliger::optimizeALS
#' @param thresh see rliger::optimizeALS
#' @param max.iters see rliger::optimizeALS
#' @param reduction.name Reduction name
#' @param reduction.key Reduction key
#' @param nrep see rliger::optimizeALS
#' @param H.init see rliger::optimizeALS
#' @param W.init see rliger::optimizeALS
#' @param V.init see rliger::optimizeALS
#' @param rand.seed Random seed
#' @param print.obj see rliger::optimizeALS
#' @param ... Arguments passed to other methods
#'
#' @importFrom Seurat DefaultAssay SplitObject GetAssayData VariableFeatures CreateDimReducObject Tool<- LogSeuratCommand
#' @importFrom rlang %||%
#'
RunOptimizeALS <- function(
  object,
  k,
  assay = NULL,
  split.by = 'orig.ident',
  lambda = 5,
  thresh = 1e-6,
  max.iters = 30,
  reduction.name = 'iNMF_raw',
  reduction.key = 'riNMF_',
  nrep = 1,
  H.init = NULL,
  W.init = NULL,
  V.init = NULL,
  rand.seed = 1,
  print.obj = FALSE,
  ...
) {

  assay <- assay %||% DefaultAssay(object = object)
  if (is.character(x = split.by) && length(x = split.by) == 1) {
    split.by <- object[[split.by]]
  }
  split.cells <- split(x = colnames(x = object), f = split.by)
  scale.data <- lapply(
    X = split.cells,
    FUN = function(x) {
      return(t(x = GetAssayData(
        object = object,
        slot = 'scale.data',
        assay = assay
      )[, x]))
    }
  )
  # scale.data <- sapply(X = scale.data, FUN = t, simplify = FALSE)
  out <- rliger::optimizeALS(
    object = scale.data,
    k = k,
    lambda = lambda,
    thresh = thresh,
    max.iters = max.iters,
    nrep = nrep,
    H.init = H.init,
    W.init = W.init,
    V.init = V.init,
    rand.seed = rand.seed,
    print.obj = print.obj
  )
  colnames(x = out$W) <- VariableFeatures(object = object)
  object[[reduction.name]] <- CreateDimReducObject(
    embeddings = do.call(what = 'rbind', args = out$H),
    loadings = t(x = out$W),
    assay = assay,
    key = reduction.key
  )
  Tool(object = object) <- sapply(
    X = out$V,
    FUN = function(x) {
      colnames(x = x) <- VariableFeatures(object = object)
      rownames(x = x) <- colnames(x = object[[reduction.name]])
      return(t(x = x))
    },
    simplify = FALSE
  )
  object <- LogSeuratCommand(object = object)
  return(object)
}


#' Run quantile_norm on a Seurat object
#'
#' @param object A merged Seurat object
#' @param split.by Attribute for splitting, defaults to "orig.ident"
#' @param reduction Embedding name
#' @param reduction.name Reduction name
#' @param reduction.key Reduction key
#' @param quantiles see rliger::quantile_norm
#' @param ref_dataset see rliger::quantile_norm
#' @param min_cells see rliger::quantile_norm
#' @param knn_k see rliger::quantile_norm
#' @param dims.use see rliger::quantile_norm
#' @param do.center see rliger::quantile_norm
#' @param max_sample see rliger::quantile_norm
#' @param eps see rliger::quantile_norm
#' @param refine.knn see rliger::quantile_norm
#' @param ... Arguments passed to other methods
#'
#' @importFrom Seurat Tool SplitObject Embeddings CreateDimReducObject DefaultAssay Tool<- Idents<- LogSeuratCommand
#'
#' @export
#'
RunQuantileNorm <- function(
  object,
  split.by = 'orig.ident',
  reduction = 'iNMF_raw',
  reduction.name = 'iNMF',
  reduction.key = 'iNMF_',
  quantiles = 50,
  ref_dataset = NULL,
  min_cells = 20,
  knn_k = 20,
  dims.use = NULL,
  do.center = FALSE,
  max_sample = 1000,
  eps = 0.9,
  refine.knn = TRUE,
  ...
) {

  embeddings <- sapply(
    X = SplitObject(object = object, split.by = split.by),
    FUN = function(x) {
      return(Embeddings(object = x[[reduction]]))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  if (is.null(x = ref_dataset)) {
    num.samples <- vapply(
      X = embeddings,
      FUN = nrow,
      FUN.VALUE = integer(length = 1L)
    )
    ref_dataset <- names(x = embeddings)[which.max(x = num.samples)]
  } else if (is.numeric(x = ref_dataset)) {
    ref_dataset <- names(x = embeddings)[ref_dataset]
  }
  if (is.character(x = ref_dataset) && !ref_dataset %in% names(x = embeddings)) {
    stop("Cannot find reference dataset '", ref_dataset, "' in the split", call. = FALSE)
  }
  out <- rliger::quantile_norm(
    object = embeddings,
    quantiles = quantiles,
    ref_dataset = ref_dataset,
    min_cells = min_cells,
    knn_k = knn_k,
    dims.use = dims.use,
    do.center = do.center,
    max_sample = max_sample,
    eps = eps,
    refine.knn = refine.knn,
    ...
  )
  object[[reduction.name]] <- CreateDimReducObject(
    embeddings = out$H.norm,
    assay = DefaultAssay(object = object[[reduction]]),
    key = reduction.key
  )
  out <- as.data.frame(x = out[names(x = out) != 'H.norm'])
  object[[colnames(x = out)]] <- out
  Idents(object = object) <- 'clusters'
  object <- LogSeuratCommand(object = object)
  return(object)
}

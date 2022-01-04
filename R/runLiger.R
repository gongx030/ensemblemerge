#' Run Seurat V3 functions
#'
#' @import Seurat
#'
#' @param data SingleCellExpirement object containing single cell counts matrix
#' @param params LigerParams object
#'
#' @return returns a SingleCellExperiment object of the integrated data
#' @export
run_Liger <- function(params, data){
  checkParams('seurat-wrappers', '0.0.0', 'R')
  
  if(class(data) != "Seurat"){
    data <- Seurat::as.Seurat(data, counts = "counts", data = NULL)
  }

  if(params@norm_data){
    print("Normalizing data")
    data <- NormalizeData(data, normalization.method = params@norm_method, 
                        scale.factor = params@scale_factor)
  }
  if(params@hvg){
    print("finding variable features")
    data <- Seurat::FindVariableFeatures(data, 
                                      selection.method = "vst", nfeatures = params@numHVG, 
                                      verbose = FALSE)
  }
  else{
      data@assays$RNA@var.features <- rownames(data)
  }
  if(params@regressUMI && params@scaling) {
    data <- ScaleData(object = data, vars.to.regress = params@vars_to_regress, do.center = FALSE, split.by = params@batch)  # in case of read count data
  } else if(params@scaling) { # default option
    data <- ScaleData(object = data, do.center = FALSE, split.by = params@batch)
  }

  data <- RunOptimizeALS(data, k = params@k, lambda = params@lambda, split.by = params@batch)
  data <- RunQuantileNorm(data, split.by = params@batch)

  if(params@return == "Seurat"){
    return(data)
  }
  else if(params@return == "SingleCellExperiment"){
    data = Seurat::as.SingleCellExperiment(data)
    return(data)
  }
  else{
    stop("Invalid return type, check params@return")
  }
}

#' Run optimizeALS on a Seurat object
#'
#' @inheritParams rliger::optimizeALS
#' @param object A merged Seurat object
#' @param split.by Attribute for splitting, defaults to "orig.ident"
#' @param ... Arguments passed to other methods
#'
#' @return A Seurat object with embeddings and loadings from \code{\link[liger]{optimizeALS}}
#' stored as a DimReduc object with name \code{reduction.name} (key set to \code{reduction.key});
#' per-dataset feature loadings matrices stored in the \code{tool} slot, accessible with
#' \code{\link[Seurat]{Tool}}
#'
#' @importFrom rliger optimizeALS
#' @importFrom Seurat DefaultAssay SplitObject GetAssayData VariableFeatures
#' CreateDimReducObject Tool<- LogSeuratCommand
#'
#' @aliases optimizeALS
#' @seealso \code{\link[liger]{optimizeALS}} \code{\link[Seurat]{Tool}}
#'
#' @export
# @method optimizeALS Seurat
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
#' @inheritParams RunOptimizeALS
#' @inheritParams rliger::quantile_norm
#' @param ... Arguments passed to other methods
#'
#' @return A Seurat object with embeddings from \code{\link[liger]{quantile_norm}}
#' stored as a DimReduc object with name \code{reduction.name} (key set to \code{reduction.key})
#'
#' @importFrom rliger quantile_norm
#' @importFrom Seurat Tool SplitObject Embeddings CreateDimReducObject
#' DefaultAssay Tool<- Idents<- LogSeuratCommand
#'
#' @aliases quantile_norm
#' @seealso \code{\link[rliger]{quantile_norm}}
#'
#' @export
# @method quantile_norm Seurat
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

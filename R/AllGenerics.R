#' Generic merge function
#' @export
setGeneric("Merge", function(x, split.by = NULL, batch_label = NULL, cell_label = NULL, merge.which = NULL, method = "Seurat", ...) standardGeneric("Merge"))
setGeneric("Score", function(x, batch_label = NULL, method = "kBET", ...) standardGeneric("Score"))
#' The generic function of Merge
#'
#' @param params a parameter object
#' @param data a data object
#' @param ... Other arguments
#'
setGeneric("Merge", function(params, data, ...) standardGeneric("Merge"))

#' The generic function of getNeighborGraph
#'
#' @param params a parameter object
#' @param data a data object
#' @param ... Other arguments
#'
setGeneric("getNeighborGraph", function(params, data, ...) standardGeneric("getNeighborGraph"))

#' The generic function of Preprocess
#'
#' @param params a parameter object
#' @param data a data object
#' @param ... Other arguments
#'
setGeneric("Preprocess", function(params, data, ...) standardGeneric("Preprocess"))

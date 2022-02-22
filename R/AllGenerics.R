#' The generic function of Merge
#'
#' @param data a data object
#' @param params a parameter object
#' @param ... Other arguments
#'
setGeneric("Merge", function(data, params, ...) standardGeneric("Merge"))

#' The generic function of getNeighborGraph
#'
#' @param data a data object
#' @param params a parameter object
#' @param ... Other arguments
#'
setGeneric("getNeighborGraph", function(data, params, ...) standardGeneric("getNeighborGraph"))

#' The generic function of Preprocess
#'
#' @param data a data object
#' @param params a parameter object
#' @param ... Other arguments
#'
setGeneric("Preprocess", function(data, params, ...) standardGeneric("Preprocess"))

#' The generic function of valid
#'
#' @param x a data object
#' @param params a parameter object
#' @param ... Other arguments
#'
setGeneric("valid", function(x, params, ...) standardGeneric("valid"))

#' The generic function of ensemble
#'
#' @param x a data object
#' @param ... Other arguments
#'
setGeneric("ensemble", function(x, ...) standardGeneric("ensemble"))

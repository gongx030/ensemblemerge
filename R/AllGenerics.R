#' The generic function of Merge
#'
#' @param data a data object
#' @param params a parameter object
#' @param ... Other arguments
#'
setGeneric("Merge", function(data, params, ...) standardGeneric("Merge"))

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
#' @param params a parameter object
#' @param ... Other arguments
#'
setGeneric("ensemble", function(x, params, ...) standardGeneric("ensemble"))

#' The generic function of Annotate
#'
#' @param x a data object
#' @param params a parameter object
#' @param ... Other arguments
#'
setGeneric("Annotate", function(x, params, ...) standardGeneric("Annotate"))

#' The generic function of Embed
#'
#' @param x a data object
#' @param params a parameter object
#' @param ... Other arguments
#'
setGeneric("Embed", function(x, params, ...) standardGeneric("Embed"))

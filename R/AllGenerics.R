#' The generic function of Merge
#'
#' @param x a data object
#' @param params a parameter object
#' @param ... Other arguments
#'
setGeneric("Merge", function(x, params, ...) standardGeneric("Merge"))

#' The generic function of Preprocess
#'
#' @param x a data object
#' @param params a parameter object
#' @param ... Other arguments
#'
setGeneric("Preprocess", function(x, params, ...) standardGeneric("Preprocess"))

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

#' The generic function of Cluster
#'
#' @param x a data object
#' @param params a parameter object
#' @param ... Other arguments
#'
setGeneric("Cluster", function(x, params, ...) standardGeneric("Cluster"))

#' The generic function of DetectDoublet 
#'
#' @param x a data object
#' @param params a parameter object
#' @param ... Other arguments
#'
setGeneric("DetectDoublet", function(x, params, ...) standardGeneric("DetectDoublet"))

#' The generic function of Normalize
#'
#' @param x a data object
#' @param params a parameter object
#' @param ... Other arguments
#'
setGeneric("Normalize", function(x, params, ...) standardGeneric("Normalize"))

#' The generic function of RemoveAmbientRNA 
#'
#' @param x a data object
#' @param params a parameter object
#' @param ... Other arguments
#'
setGeneric("RemoveAmbientRNA", function(x, params, ...) standardGeneric("RemoveAmbientRNA"))

#' The generic function of ReferenceMap 
#'
#' @param query a data object
#' @param atlas a data object
#' @param params a parameter object
#' @param ... Other arguments
#'
setGeneric("ReferenceMap", function(query, atlas, params, ...) standardGeneric("ReferenceMap"))


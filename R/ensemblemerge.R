#' Get valid integration methods 
#'
#'
#' @return returns vector of valid names for integration methods
#' @export
getMethods <- function(){
  Methods = c("Seurat", "Harmony", "Liger", "Scanorama", "BBKNN", "fastMNN", "Uncorrected")
  return(Methods)
}

#' Get valid parameter classes 
#'
#'
#' @return returns vector of valid names for integration methods
getParams <- function(){
  Params = c("SeuratParams", "HarmonyParams", "LigerParams", "ScanoramaParams", "BBKNNParams", "FastMNNParams", "UncorrectedParams")
  return(Params)
}


#' Run enhsemblemerge functions
#'
#' @import Seurat
#' @import SingleCellExperiment
#'
#' @param x SummarizedExpirement object containing single cell counts matrix
#' @return returns a SummarizedExperiment object of the integrated data
#' @export
EnsembleMerge <- function(data, methods = c("Seurat", "Harmony"), return = "SingleCellExperiment"){

  #* Add error if methods is not of class character or xParams
  if(class(methods) != "character" & class(methods) != "list"){
      stop("'methods' must be of class character or list of of Params objects set by setParams()")
  }

  #* Add error if methods length is not >= 2
  if(length(methods) < 2){
    stop("'methods' must be a list or vector of at least 2 valid integration methods")
  }

  if(class(methods[1][[1]]) %in% getParams()){
    #* Add error if invalid methods present
    for(i in 1:length(methods)){
      if(!(methods[i][[1]] %in% getParams())){
        stop(sprintf("'%s' is not an available method, check available methods with getParams()", methods[i]))
      }
    }
  }
  else{
    new_methods = list() #convert methods to params objects
    #* Add error if invalid methods present, also convert methods to params objects if no error
    for(i in 1:length(methods)){
      if(!(methods[i] %in% getMethods())){
        stop(sprintf("'%s' is not an available method, check available methods with getMethods()", methods[i]))
      }
      else{
        new_methods[[i]] = setParams(method = methods[i]) #assign params object to methods
      }
    }
    methods = new_methods #exchange variable name for downstream
    new_methods = NULL #remove temp methods holder
  }
  

  #* Add error if data is not a valid input type
  if(!(class(data) %in% c("Seurat", "SummarizedExperiment", "SingleCellExperiment"))){
      stop(sprintf("'%s' is not an valid input for data, available inputs are 'Seurat', 'SummarizedExperiment', 'SingleCellExperiment'", class(data)))
  }
  
  if(class(data) == 'SummarizedExperiment'){
      data = as(data, "SingleCellExperiment") # set data as SingleCellExperiment
    }

    if(class(data) == "Seurat"){
      data = Seurat::as.SingleCellExperiment(data)
    }

  ng = lapply(methods, function(ng){
    ng = getNeighborGraph(ng, data) 
  })

  ng = Reduce("+", ng)/length(ng)

  if(return == "SingleCellExperiment"){
    S4Vectors::metadata(data)$EnsembleMerge <- ng
  }
  else if(return == "Seurat"){
    data = as.Seurat(data, counts = "counts", data = NULL)
    data[["EnsembleMerge"]] = as.Graph(ng)
  }
  

  return(data)
}

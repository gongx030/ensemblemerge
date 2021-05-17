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


#' Run ensemblemerge functions
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

  agreement = function(x){
    #ret = Matrix::Matrix(nrow = nrow(x[[1]]), ncol = ncol(x[[1]]), data = 0, sparse = TRUE)
    #ret = as(ret, "dgCMatrix")
    #ret = as(ret, "dgCMatrix")
    ret = ng[[1]]
    temp = Reduce(rbind, lapply(x, summary))
    temp = temp[which(duplicated(temp[,c("i", "j")])),]
    temp = unique(temp[,c("i", "j")])
    message("applying boolean values")
    ret[as.matrix(temp[,c("i", "j")])] = 1
    return(ret)
  }

  #sigmoid = suppressWarnings(function(x){0.5 * (1 + sin((x*pi)-(pi/2)))})
  sigmoid = function(x){1/(1+(x/(1-x))^-5)}

  agreement = agreement(ng)

  wt = list()
  kg = list()
  for(i in 1:length(ng)){
    weight = sum(ng[[i]]*agreement)/sqrt((ng[[i]]*ng[[i]])%*%(agreement*agreement))
    x = ng[[i]]*weight
    x[is.na(x)] = 0
    weight[is.na(weight)] = 0
    x = as(x, "dgCMatrix")
    weight = as(weight, "dgCMatrix")
    wt = append(wt, weight)
    kg = append(kg, x)
  }

tmpwt = as.vector(unlist(wt))
for(i in 1:length(wt)){
  wt[i] = (wt[[i]] - min(tmpwt))/(max(tmpwt) - min(tmpwt))
  message(sprintf("weight is %f", wt[[i]]))
  wt[i] = sigmoid(wt[[i]])
  message(sprintf("weight is now %f", wt[i]))
}

for(i in 1:length(methods)){
    message(sprintf("%s method scores %f out of 1", methods[[i]]@name, wt[[i]]))
}

  ng = Reduce("+", kg)/Reduce("+", wt)
  ng[is.na(ng)] = 0
  ng = as(ng, "dgCMatrix")

  if(return == "SingleCellExperiment"){
    S4Vectors::metadata(data)$EnsembleMerge <- ng
  }
  else if(return == "Seurat"){
    data = as.Seurat(data, counts = "counts", data = NULL)
    data[["EnsembleMerge"]] = as.Graph(ng)
  }
  

  return(data)
}

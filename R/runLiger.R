#' Run LIGER functions
#'
#' @import rliger
#' @import Seurat
#'
#' @param x SummarizedExpirement object containing single cell counts matrix
#' @return returns a SummarizedExperiment object of the integrated data
#' @export
liger_preprocess <- function(x, batch_label = "batchlb", var.thresh=0.1)
{
  metadata = colData(x)
  expr_mat = assay(x)
  data <- lapply(levels(metadata[,batch_label]), function(l){
    cells <- metadata$cell[metadata[,batch_label]==l]
    return(expr_mat[,cells])
  })
  names(data) <- levels(metadata[,batch_label])
  ligerex <- createLiger(data)
  
 
  # Normalize the data to account for different numbers of UMIs per cell 
  ligerex = normalize(ligerex)  
  
  # Select variable genes on each of the datasets separately, then takes the union
  ligerex = selectGenes(ligerex, var.thresh = var.thresh)
  
  # Scale the data
  # note : as.matrix() does not work if matrix has more than 200000 rows  
  
  dim_data <- sapply(data,function(l){dim(l)[1]})
  
  chunk <- function(x,n){
    vect <- c(1:x)
    num <- ceiling(x/n)
    split(vect,rep(1:num,each=n,len=x))
  }
  
  if(any(dim_data>200000)){
    ligerex@scale.data <- lapply(1:length(ligerex@norm.data),function(i) {liger:::scaleNotCenterFast(t(ligerex@norm.data[[i]][ligerex@var.genes,]))})
    ligerex@scale.data <- lapply(ligerex@scale.data,function(l){
      if(dim(l)[1]>200000){
        l2 <- lapply(chunk(nrow(l),200000), function(i){as.matrix(l[i,])})
        res <- do.call(rbind,l2)
      } else {
        res <- as.matrix(l)
      }
      return(res)
    })
    names(ligerex@scale.data) <- names(ligerex@norm.data)
    for (i in 1:length(ligerex@scale.data)) {
      ligerex@scale.data[[i]][is.na(ligerex@scale.data[[i]])] <- 0
      rownames(ligerex@scale.data[[i]]) <- colnames(ligerex@raw.data[[i]])
      colnames(ligerex@scale.data[[i]]) <- ligerex@var.genes
    }
    ligerex <- liger:::removeMissingObs(ligerex, slot.use = "scale.data", use.cols = F)
  } else {
    ligerex = scaleNotCenter(ligerex)
  }
  
  return(ligerex)
}

#' @export
call_liger <- function(ligerex, batch_label,
                       k = 20, nrep = 3)
{
  
  ligerex = optimizeALS(ligerex, k = k, lambda = 5, nrep = nrep)
    
  ligerex = quantileAlignSNF(ligerex)

  return(ligerex)
}

#' @export
run_LIGER <- function(x, batch_label = "batchlb", var.thresh=0.1, k = 20, nrep = 3){
  ligerex <- liger_preprocess(expr_mat, batch_label, var.thresh=0.1)
  integrated <- call_liger <- function(ligerex, batch_label, k, nrep)
  integrated <- ligerToSeurat(integrated)
  integrated <- as.SingleCellExperiment(integrated)
  return(integrated)
}
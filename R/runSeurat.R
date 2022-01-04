#' Run Seurat V3 functions
#'
#' @import Seurat
#'
#' @param x SummarizedExpirement object containing single cell counts matrix
#' @return returns a SummarizedExperiment object of the integrated data
#' @export
seurat3_preprocess <- function(x, 
                              normData = TRUE, Datascaling = TRUE, 
                              min_cells = 10, min_genes = 300, 
                              norm_method = "LogNormalize", 
                              scale_factor = 10000, 
                              numVG = 300, nhvg = 2000, 
                              batch_label = "batchlb", celltype_label = "CellType",
                              hvg = T, k.weight = 100)
{
  checkParams('Seurat', '4.0.1', 'R')
  ##########################################################
  # preprocessing
  batches = x
  if(class(batches) != "Seurat"){
    batches <- Seurat::as.Seurat(batches, counts = "counts", data = NULL)
  }

  batch_list <- Seurat::SplitObject(batches, split.by = batch_label)
  for(i in 1:length(batch_list)) {
    if(normData == TRUE){
      batch_list[[i]] <- NormalizeData(batch_list[[i]],normalization.method = norm_method, scale.factor = scale_factor)
    }
    if(hvg){
      batch_list[[i]] <- Seurat::FindVariableFeatures(batch_list[[i]], selection.method = "vst", nfeatures = nhvg, 
                                              verbose = FALSE)
      #print(dim(batch_list[[i]]))
    }
    else{
      batch_list[[i]]@assays$RNA@var.features <- rownames(batch_list[[i]])
    }
  }
  return(batch_list)
}
#' @export
call_seurat3 <- function(batch_list, batch_label, celltype_label, npcs = 20, seed = 1, regressUMI = TRUE, Datascaling = TRUE, dims = 20, pca_name = "PCA", k.weight = 100)
{
  for(batch in batch_list){
    if(nrow(batch) < npcs){
      stop(sprintf("batch: %s has %i cells, all batches must have more cells than dimensions for anchor weighting (dims:%i)", unique(batch@meta.data[[batch_label]]), nrow(batch), dims))
    }
  }

  #check if batches are less than k.weight
  for(batch in batch_list){
    if(ncol(batch) < k.weight){
      stop(sprintf("batch: %s has %i cells, all batches must have more cells than the value of k.weight (k.weight:%i)", unique(batch@meta.data[[batch_label]]), ncol(batch), k.weight))
    }
  }

  cell_anchors <- Seurat::FindIntegrationAnchors(object.list = batch_list, dims = 1:dims)
  batches <- Seurat::IntegrateData(anchorset = cell_anchors, dims = 1:npcs, k.weight = k.weight)
  dim(batches)

  Seurat::DefaultAssay(batches) <- "integrated"

  if(regressUMI && Datascaling) {
    batches <- Seurat::ScaleData(object = batches, vars.to.regress = c("nUMI"))  # in case of read count data
  } else if(Datascaling) { # default option
    batches <- Seurat::ScaleData(object = batches)
  }

  batches <- Seurat::RunPCA(object = batches, npcs = npcs, verbose = FALSE, reduction.name = pca_name)

  return(batches)
}
#' @export
run_Seurat <- function(params, data){
                        batch_list = seurat3_preprocess(x = data, 
                              normData = params@norm_data, Datascaling = params@scaling, 
                              min_cells = params@min_cells, min_genes = params@min_genes, 
                              norm_method = params@norm_method, 
                              scale_factor = params@scale_factor, 
                              numVG = params@numVG, nhvg = params@numHVG, 
                              batch_label = params@batch, k.weight = params@k.weight)
                        integrated = call_seurat3(batch_list, batch_label = params@batch, 
                                                  npcs = params@npcs, seed = params@seed, 
                                                  regressUMI = params@regressUMI, Datascaling = params@scaling,
                                                  dims = params@dims, pca_name = params@dimreduc_names[["PCA"]],
                                                  k.weight = params@k.weight)

                        
                        if(params@return == "Seurat"){
                          return(integrated)
                        }
                        else if(params@return == "SingleCellExperiment"){
                          RNA = Seurat::as.SingleCellExperiment(integrated, assay = "RNA")
                          integrated = Seurat::as.SingleCellExperiment(integrated)
                          SingleCellExperiment::altExp(integrated, params@altExp_names) = RNA
                          return(integrated)
                        }
                        else{
                          stop("Invalid return type, check params@return")
                        }
                      }

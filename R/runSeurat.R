#' Run Seurat V3 functions
#'
#' @import Seurat
#'
#' @param x SummarizedExpirement object containing single cell counts matrix
#' @return returns a SummarizedExperiment object of the integrated data
#' @export
seurat3_preprocess <- function(x, 
                              normData = TRUE, Datascaling = TRUE, regressUMI = FALSE, 
                              min_cells = 10, min_genes = 300, 
                              norm_method = "LogNormalize", 
                              scale_factor = 10000, 
                              numVG = 300, nhvg = 2000, 
                              batch_label = "batchlb", celltype_label = "CellType",
                              hvg = T)
{

  ##########################################################
  # preprocessing

  batches <- Seurat::as.Seurat(x, counts = "counts", data = NULL)

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
      batch_list[[i]]@assays$RNA@var.features <- rownames(expr_mat_seurat)
    }
  }
  return(batch_list)
}
#' @export
call_seurat3 <- function(batch_list, batch_label, celltype_label, npcs = 20, seed = 1, regressUMI = TRUE, Datascaling = TRUE)
{
  cell_anchors <- Seurat::FindIntegrationAnchors(object.list = batch_list, dims = 1:npcs)
  batches <- Seurat::IntegrateData(anchorset = cell_anchors, dims = 1:npcs)
  dim(batches)

  Seurat::DefaultAssay(batches) <- "integrated"

  if(regressUMI && Datascaling) {
    batches <- Seurat::ScaleData(object = batches, vars.to.regress = c("nUMI"))  # in case of read count data
  } else if (Datascaling) { # default option
    batches <- Seurat::ScaleData(object = batches)
  }

  batches <- Seurat::RunPCA(object = batches, npcs = npcs, verbose = FALSE)

  batches <- Seurat::RunUMAP(batches, reduction = "pca", dims = 1:npcs, k.seed = seed)
}
#' @export
run_Seurat <- function(x, 
                      normData = TRUE, Datascaling = TRUE, regressUMI = FALSE, 
                      min_cells = 10, min_genes = 300, 
                      norm_method = "LogNormalize", 
                      scale_factor = 10000, 
                      numVG = 300, nhvg = 2000, 
                      batch_label = "batchlb", celltype_label = "CellType",
                      hvg = T, npcs = 20, seed = 1){
                        batch_list = seurat3_preprocess(x, 
                              normData, Datascaling, regressUMI, 
                              min_cells, min_genes, 
                              norm_method, 
                              scale_factor, 
                              numVG, nhvg, 
                              batch_label, celltype_label)
                        integrated = call_seurat3(batch_list, batch_label, celltype_label, npcs, seed, regressUMI, Datascaling)
                        integrated = Seurat::as.SingleCellExperiment(integrated)
                        return(integrated)
                      }

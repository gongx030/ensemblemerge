#' Run Harmony functions
#'
#' @import harmony
#' @import Seurat
#'
#' @param x SingleCellExperiment object containing single cell counts matrix
#' @return returns a SingleCellExperiment object of the integrated data
#' @export
harmony_preprocess <- function(x,
                              normData = TRUE, Datascaling = TRUE, regressUMI = FALSE, 
                              norm_method = "LogNormalize", scale_factor = 10000, 
                              nfeatures = 300, npcs = 20, 
                              batch_label = "batchlb", celltype_label = "CellType")
{

  b_seurat <- as.Seurat(x, counts = "counts", data = NULL)

  if(normData){
    b_seurat <- NormalizeData(object = b_seurat, normalization.method = norm_method, scale.factor = scale_factor)
  } else{
    b_seurat@data = b_seurat@raw.data
  }

  b_seurat <- Seurat::FindVariableFeatures(object = b_seurat, nfeatures = nfeatures)

  if(regressUMI && Datascaling) {
    b_seurat <- ScaleData(object = b_seurat, vars.to.regress = c("nUMI"))  # in case of read count data
  } else if (Datascaling) { # default option
    b_seurat <- ScaleData(object = b_seurat)
  }

  return(b_seurat)
}

#' @export
call_harmony_2 <- function(b_seurat, batch_label, celltype_label, npcs = 20, seed = 1)
{

  #Harmony settings
  theta_harmony = 2
  numclust = 50
  max_iter_cluster = 100


  b_seurat <- RunPCA(object = b_seurat, npcs = npcs, pc.genes = b_seurat@var.genes, verbose = FALSE)

  b_seurat <- RunHarmony(object = b_seurat, batch_label, theta = theta_harmony, plot_convergence = TRUE, 
                          nclust = numclust, max.iter.cluster = max_iter_cluster)


  b_seurat <- RunUMAP(b_seurat, reduction.use = "harmony", dims = 1:npcs, k.seed = seed)

  return(b_seurat)
}

#' @export
run_Harmony <- function(x,
                        normData = TRUE, Datascaling = TRUE, regressUMI = FALSE, 
                        norm_method = "LogNormalize", scale_factor = 10000, 
                        numVG = 300, npcs = 20, 
                        batch_label = "batchlb", celltype_label = "CellType",
                        seed = 1){
  b_seurat = harmony_preprocess(x,
                              normData = TRUE, Datascaling = TRUE, regressUMI = FALSE, 
                              norm_method, scale_factor, 
                              numVG, npcs, 
                              batch_label, celltype_label)
  integrated = call_harmony_2(b_seurat, batch_label, celltype_label, npcs, seed)
  integrated = as.SingleCellExperiment(integrated)
  return(integrated)
}

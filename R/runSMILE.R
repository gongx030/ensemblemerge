#' Run Seurat V3 functions
#'
#' @import reticulate
#' @import Seurat
#'
#' @param data SingleCellExperiment object containing single cell counts matrix
#' @param params SMILEParams object generated from setParams(method = "BBKNN") function
#' @return returns a SummarizedExperiment object of the integrated data
#' @export
run_SMILE <- function(params, data){
  ### load up python environment ###
  reticulate::py_config()
  sc <<- reticulate::import("scanpy", delay_load = TRUE)
  sr <<- reticulate::import("scanorama", delay_load = TRUE)
  ad <<- reticulate::import("anndata", delay_load = TRUE, convert = FALSE)
  bbknn <<- reticulate::import("bbknn", delay_load = TRUE)
  #filepath = system.file("R/runSMILE.py", package = "ensemblemerge")
  py$adata = suppressWarnings(sceasy::convertFormat(data, from = "sce", to = "anndata"))
  py_run_string("import numpy as np
import pandas as pd
import scanpy as sc
from SMILE import SMILE
from SMILE.SMILE import SMILE_trainer

import torch

import umap
import anndata
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics import normalized_mutual_info_score

import warnings
warnings.filterwarnings('ignore')

import seaborn as sns

def weights_init(m):
    if isinstance(m, torch.nn.Linear):
        torch.nn.init.xavier_uniform(m.weight.data)
        m.bias.data.zero_()


X = adata.X
scaler = StandardScaler()
X = scaler.fit_transform(X)

cells = adata.obs['cell'].values.tolist()
cells = np.array(cells)


X_all_tensor = torch.tensor(X).float()
clf_out = 25
net = SMILE.SMILE(input_dim=X.shape[1],clf_out=clf_out)
net.apply(weights_init)

SMILE_trainer(X,net)

net.to(torch.device('cpu'))
y_pred = net.encoder(X_all_tensor)
y_pred = torch.Tensor.cpu(y_pred).detach().numpy()

embedding = umap.UMAP(n_neighbors=10, min_dist=0.25, n_components=2,
                        metric='euclidean').fit_transform(y_pred)
embedding = pd.DataFrame(embedding)
embedding.columns=['UMAP1','UMAP2']
embedding['Cluster']=cells
adata.obsm['X_umap']=embedding.iloc[:,:2].values
adata.obs['CellType']=cells
adata.write(filename = 'temp.h5ad')")
  #source_python(filepath)
  integrated = sceasy::convertFormat("temp.h5ad", from = "anndata", to = "seurat")
  return(integrated)
}
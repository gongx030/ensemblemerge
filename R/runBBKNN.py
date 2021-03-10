import numpy as np
import bbknn
import scanpy as sc

sc.pp.filter_cells(adata, min_genes=300)
sc.pp.filter_genes(adata, min_cells=5)

sc.pp.log1p(adata)
sc.pp.scale(adata)
sc.tl.pca(adata, svd_solver=svd_solver)
sc.pp.neighbors(adata,n_neighbors=int(n_neighbors), n_pcs=int(npcs))

adata_bbknn = bbknn.bbknn(adata, copy=copy, neighbors_within_batch=int(neighbors_within_batch),
                          approx=approx,trim=int(trim))

sc.tl.pca(adata_bbknn, svd_solver=svd_solver,n_comps=int(npcs))
adata_bbknn.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat
sc.tl.umap(adata_bbknn)
adata_bbknn.write(filename = "temp.h5ad")
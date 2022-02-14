<img src = "https://github.com/gongx030/gongx030.github.io/blob/cac14d855edd4e785056c09e340e5213c7b6e572/images/ensemblemerge_graphics_v2.png?raw=true" width = "250" height = "200" align = "right" />

# ensemblemerge

<!-- badges: start -->
[![R-CMD-check](https://github.com/erikjskie/ensemblemerge/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/erikjskie/ensemblemerge/actions/workflows/check-standard.yaml)
[![Github tag](https://badgen.net/github/tag/Naereen/Strapdown.js)](https://github.com/erikjskie/ensemblemerge/tags/)
<!-- badges: end -->


ensemblemerge is a package that implements a common work flow for several single cell RNA sequence integration methods including Seurat, Scanorama, Harmony, Liger, fastMNN, bbknn and scVI. The merging process is designed to be as streamlined as possible for the user, only requiring a **SingleCellExperiment**, or **Seurat** object and a specification of which method should be used to integrate in `setParams()`.

## Prerequisites

### Step 1. We recommend to build a new `conda` environment for ensemblemerge and install several R packages required for base function:

```
conda create -n ensemblemerge python=3.7
conda activate ensemblemerge
```

```
conda install -y -c conda-forge r-base=4.1.2
conda install -y -c conda-forge r-devtools=2.4.2
conda install -y -c conda-forge r-seurat=4.1.0
conda install -y -c conda-forge r-r.utils=2.11.0
conda install -y -c conda-forge umap-learn=0.5.2
```

### Step 2: Installing ensemblemerge:

```r
devtools::install_github("erikjskie/ensemblemerge")
```

### Step 3: Installing dependencies:

To use the constituent methods other than `Seurat`, the users will need to have the R or Python dependencies installed manually. For R packages, ensemblemerge uses `packageVersion` to verify required packages and versions. For Python packages, ensemblemerge ueses the system command `pip show` to verify required packages and versions. 

| Method | R | Python |
| --- | --- | --- |
| Seurat | | |
| Harmony | `install.packages("harmony")` | |
| LIGER | `install.packages('rliger')` | |
| fastMNN | `BiocManager::install('batchelor')` | |
| Scanorama |  | `pip install scanpy==1.8.2` <br> `pip install anndata==0.7.8` <br> `pip install scanorama==1.7.1`  | 
| scVI |  | `pip install scanpy==1.8.2` <br> `pip install anndata==0.7.8` <br> `pip install scvi-tools==0.14.5` | 
| BBKNN |  | `pip install bbknn==1.5.1` <br> `pip install scanpy==1.8.2` <br> `pip install anndata==0.7.8` <br> `pip install leidenalg==0.8.8 ` |


## Quick start guide

Once ensemblemerge is installed, merging batches can be performed by the following:

```r
library(ensemblemerge)

#merge with single method
pbmc #small 500 cell by 500 feature, 2 batch dataset

params = setParams(method = "Seurat", return = "Seurat") #different integration methods can be selected by setting method, see methods by calling getMethods()

merged_data = Merge(params, pbmc)

ensemblemerged_data = EnsembleMerge(pbmc, methods = c("Seurat", "Harmony", "BBKNN"), return = "Seurat")
```

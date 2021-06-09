[![Build Status](https://travis-ci.com/erikjskie/ensemblemerge.svg?token=TzArZ5EDidamcqdAtCie&branch=main)](https://travis-ci.com/github/erikjskie/ensemblemerge)
[![codecov](https://codecov.io/gh/erikjskie/ensemblemerge/branch/main/graph/badge.svg?token=J1H0OAEQ5S)](https://codecov.io/gh/erikjskie/ensemblemerge)
# ensemblemerge
ensemblemerge is a package that implements a common work flow for several single cell RNA sequence integration methods including Seurat V3, Scanorama, Harmony, Liger, fastMNN, and bbknn. The merging process is designed to be as streamlined as possible for the user, only requiring a **SummarizedExperiment**, **SingleCellExperiment**, or **Seurat** object and a specification of which method should be used to integrate in `setParams()`.

**Instructions, documentation, and examples can be found online at: [EnsembleMerge](https://erikjskie.github.io/packages/ensemblemerge/)**

## Prerequisites
The following R packages are required for installation of ensemblemerge:

```{R}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("SingleCellExperiment", "SummarizedExperiment", "LoomExperiment"))
devtools::install_github('satijalab/seurat-wrappers')
devtools::install_github("cellgeni/sceasy")
```

The following packages are necessary only if using the following Merge methods:
  * For Harmony:
  ```{R}
  devtools::install_github("immunogenomics/harmony")
  ```
  * For Scanorama:
  ```{Python}
  pip3 install scanpy
  pip3 install scanorama
  ```
  * For bbknn:
  For R
  ```{R}
  remotes::install_github("rstudio/reticulate") #increases available memory
  ```
  Python Packages
  ```{Python}
  system("pip3 install pynndescent") #optimizes dimension reduction
  system("pip3 install leidenalg") #optional clustering algorithm that improves bbknn performance
  ```
## Quick start guide

Once ensemblemerge is installed, merging batches can be performed by the following:

```{R}
library(ensemblemerge)

data #some SingleCellExperiment or SummarizedExperiment object with batches to be integrated

params = setParams() #different integration methods can be selected by setting method = c("Seurat", "Harmony", "Liger", "Scanorama", "BBKNN")

merged_data = Merge(params, data)
```

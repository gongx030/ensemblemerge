# ensemblemerge

<!-- badges: start -->
[![R-CMD-check](https://github.com/erikjskie/ensemblemerge/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/erikjskie/ensemblemerge/actions/workflows/check-standard.yaml)
<!-- badges: end -->

## 1. A simple scRNA-seq pipeline

```
library(ensemblemerge)
x # a Seurat object

# preprocessing
params_preprocess <- new('SeuratPreprocess')
x <- Preprocess(x, params_preprocess)

# normalization
params_normalize <- new('SeuratNormalize', preprocess = params_preprocess)
x <- Normalize(x, params_normalize)

# doublet removing 
params_doubletdetect <- new('DoubletFinderDoubletDetect',  normalize = params_normalize)
x <- DetectDoublet(x, params_doubletdetect)

# dimension reduction
params_embed <- new('PCAEmbed', normalize = params_normalize)
x <- Embed(x, params_embed)

# clustering
params_cluster <- new('LouvainCluster', embedding = params_embed)
x <- Cluster(x, params_cluster)

# annotation
# params_genemarkers <- new('PanglaoDBGeneMarkers', genome = 'hg19')
# params_annotate <- new('clustifyrAnnotate', normalize = params_normalize, gene_marker = params_genemarkers, cluster = params_cluster)
# x <- Annotate(x, params_annotate)
```

## 2. An integration pipeline
```
library(ensemblemerge)
x # a Seurat object

# preprocessing
params_preprocess <- new('SeuratPreprocess', batch = 'batch')
x <- Preprocess(x, params_preprocess)

# normalization
params_normalize <- new('SeuratNormalize', preprocess = params_preprocess)
x <- Normalize(x, params_normalize) # the returned x is a SeuratList object

# doublet removing 
params_doubletdetect <- new('scDblFinderDoubletDetect',  normalize = params_normalize)
x <- DetectDoublet(x, params_doubletdetect)

# integration
params_merge <- new('SeuratMerge', normalize = params_normalize)
x_merged <- Merge(x, params_merge) # x_merged is a Seurat object
```

## Preprocessing

* cell and gene filtering 
* Check the integrity of the input data.  For example, if a `batch` indicator is specified, it will validate whether the meta data is available. 
* No batch related operations

| Class | Method | Ref. |
| --- | --- | --- |
| `SeuratPreprocess` |  |  | 



## Normalization
| Class | Method | Ref. |
| --- | --- | --- |
| `SCTransformNormalize` | [SCTransform](https://satijalab.org/seurat/articles/sctransform_vignette.html) | [[Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1)] | 
| `SeuratNormalize` | [Seurat](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)'s default normalization pipeline | [[Paper](https://www.nature.com/articles/nbt.4096)] | 

## Doublet detection

* The detected doublets are removed. 
* If a `SeuratList` object is provided, the doublet detection is performed on each batch separately. 

| Class | Method | Ref. |
| --- | --- | --- |
| `scDblFinderDoubletDetect` | [scDblFinder](https://bioconductor.org/packages/release/bioc/html/scDblFinder.html) | [[Paper](https://f1000research.com/articles/10-979)] | 
| `DoubletFinderDoubletDetect` | [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) | [[Paper](https://pubmed.ncbi.nlm.nih.gov/30954475/)] | 

## Ambient RNA removal

* The adjusted counts will be used to replace the raw counts. 

| Class | Method | Ref. |
| --- | --- | --- |
| `decontXRemoveAmbientRNA` | [celda](http://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html) | [[Yang et al. 2020](https://doi.org/10.1186/s13059-020-1950-6)] | 


## De novo integration
| Class | Method | Ref. |
| --- | --- | --- |
| `SeuratMerge` | [Seurat(v3)](https://satijalab.org/seurat/articles/integration_introduction.html) |  | 
| `HarmonyMerge` | [Harmony](https://github.com/immunogenomics/harmony) | [[Paper](https://www.nature.com/articles/s41592-019-0619-0)]
| `FastMNNMerge` | [FastMNN](https://rdrr.io/github/LTLA/batchelor/man/fastMNN.html) | [[Paper](https://www.nature.com/articles/nbt.4091)]
| `LigerMerge` | [LIGER](https://www.nature.com/articles/s41596-020-0391-8) | [[Paper](https://www.nature.com/articles/s41596-020-0391-8)]
| `BBKNNMerge` | | |
| `ScanoramaMerge` | | |
| `scVIMerge` | | |



# ensemblemerge

<!-- badges: start -->
[![R-CMD-check](https://github.com/erikjskie/ensemblemerge/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/erikjskie/ensemblemerge/actions/workflows/check-standard.yaml)
<!-- badges: end -->

## 1. A simple scRNA-seq pipeline

```
library(ensemblemerge)
x # a Seurat object

# preprocessing
params_preprocess <- new('SeuratPreprocess', batch = 'batch')
x <- Preprocess(x, params_preprocess)

# normalization
params_normalize <- new('SeuratNormalize', preprocess = params_preprocess)
x <- Normalize(x, params_normalize)

# doublet removing 
params_doubletdetect <- new('DoubletFinderoubletDetect',  normalize = params_normalize)
x <- DetectDoublet(x, params_doubletdetect)

# dimension reduction
params_embed <- new('PCAEmbed', normalize = params_normalize)
x <- Embed(x, params_doubletdetect)

# clustering
params_cluster <- new('LouvainCluster', embedding = params_embed)
x <- Cluster(x, params_cluster)

# annotation
# params_genemarkers <- new('PanglaoDBGeneMarkers', genome = 'hg19')
# params_annotate <- new('clustifyrAnnotate', normalize = params_normalize, gene_marker = params_genemarkers, cluster = params_cluster)
# x <- Annotate(x, params_annotate)
```

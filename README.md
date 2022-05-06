# ensemblemerge

<!-- badges: start -->
[![R-CMD-check](https://github.com/erikjskie/ensemblemerge/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/erikjskie/ensemblemerge/actions/workflows/check-standard.yaml)
<!-- badges: end -->

## Examples

| Task | Colab | Jupyter | Version(Tag) |
| --- | --- | --- | --- |
| scRNA-seq normalization | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/gongx030/ensemblemerge_docs/blob/main/notebooks/normalization.ipynb) | [![Jupyter Notebook](https://img.shields.io/badge/jupyter-%23FA0F00.svg?style=for-the-badge&logo=jupyter&logoColor=white)](https://github.com/gongx030/ensemblemerge_docs/blob/main/notebooks/normalization.ipynb) | `v2.1.21-001` |
| Reference-based scRNA-seq integration | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/gongx030/ensemblemerge_docs/blob/main/notebooks/reference_based_integration.ipynb) | [![Jupyter Notebook](https://img.shields.io/badge/jupyter-%23FA0F00.svg?style=for-the-badge&logo=jupyter&logoColor=white)](https://github.com/gongx030/ensemblemerge_docs/blob/main/notebooks/reference_based_integration.ipynb) | `v2.1.23-001` |
| *de novo* scRNA-seq integration | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/gongx030/ensemblemerge_docs/blob/main/notebooks/de_novo_integration_v1.ipynb) | [![Jupyter Notebook](https://img.shields.io/badge/jupyter-%23FA0F00.svg?style=for-the-badge&logo=jupyter&logoColor=white)](https://github.com/gongx030/ensemblemerge_docs/blob/main/notebooks/de_novo_integration_v1.ipynb) | `v2.1.24-001` |


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
| `ScrubletDoubletDetect` | [Scrublet](https://github.com/swolock/scrublet) | [[Wolock et al. 2019](https://pubmed.ncbi.nlm.nih.gov/30954476/)] | 

## Ambient RNA removal

* The adjusted counts will be used to replace the raw counts. 

| Class | Method | Ref. |
| --- | --- | --- |
| `decontXRemoveAmbientRNA` | [celda](http://bioconductor.org/packages/release/bioc/vignettes/celda/inst/doc/decontX.html) | [[Yang et al. 2020](https://doi.org/10.1186/s13059-020-1950-6)] | 
| `soupXRemoveAmbientRNA` | [SoupX](https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html) | [[Young et al. 2020](https://academic.oup.com/gigascience/article/9/12/giaa151/6049831?login=true)] | 

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



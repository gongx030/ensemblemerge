# syntax=docker/dockerfile:1
FROM mambaorg/micromamba:latest

RUN micromamba install --yes --name base --channel conda-forge --channel bioconda --channel r \
			python=3.7 \
			r-base=4.2 \
	   	r-seurat=4.3 \
    	scanpy>=1.9 \
    	scvi-tools>=0.19 \
    	numpy \
    	anndata \
	    r-matrix \
 	   	r-dplyr \
	    r-magrittr \
	    r-rlang \
	    r-uwot \
	    r-harmony \
	    r-liger \
	    r-sceasy \
			r-scina \
    	r-rmtstat \
    	r-mclust \
    	r-cluster \
    	r-clue \
    	r-parallelly \
#    	r-soupx \
    	bioconductor-singlecellexperiment \
    	bioconductor-summarizedexperiment \
    	bioconductor-basilisk \
    	bioconductor-biocgenerics \
    	bioconductor-zellkonverter \
    	bioconductor-sc3 \
    	bioconductor-simlr \
    	bioconductor-scran \
    	bioconductor-clustifyr \
    	bioconductor-scdblfinder \
 	  	bioconductor-celda \
	    bioconductor-scuttle \
	   	bioconductor-s4vectors \
			r-reticulate && \
		micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1  

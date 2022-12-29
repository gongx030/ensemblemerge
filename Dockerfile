# syntax=docker/dockerfile:1
FROM --platform=linux/amd64 mambaorg/micromamba:latest

RUN  apt-get update \
  && apt-get install -y wget \
  && rm -rf /var/lib/apt/lists/*

RUN micromamba install --yes --name base --channel conda-forge --channel bioconda --channel r \
			python=3.7 \
			r-base=4.2.2 \
	   	r-seurat=4.3 \
    	numpy>=1.21 \
    	scanpy>=1.9 \
    	scvi-tools>=0.19 \
    	anndata>=0.8 \
	    r-matrix \
 	   	r-dplyr \
	    r-magrittr \
	    r-rlang \
	    r-uwot \
	    r-harmony \
	    r-liger >=0.5.0 \
	    r-sceasy >=0.0.6 \
			r-scina >=1.2.0 \
    	r-rmtstat \
    	r-mclust >=5.4.9 \
    	r-cluster >=2.1.2 \
    	r-clue >=0.3.60 \
    	r-parallelly \
    	r-soupx >=1.5.2 \
    	bioconductor-singlecellexperiment \
    	bioconductor-summarizedexperiment \
    	bioconductor-basilisk \
    	bioconductor-biocgenerics \
    	bioconductor-zellkonverter \
    	bioconductor-sc3 >=1.22.0 \
    	bioconductor-simlr >=1.18.0 \
    	bioconductor-scran >=1.22.1 \
    	bioconductor-clustifyr \
    	bioconductor-scdblfinder >=1.8.0 \
 	  	bioconductor-celda >=1.10.0 \
	    bioconductor-scuttle >=1.4.0 \
	   	bioconductor-s4vectors \
			r-reticulate && \
		micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1  

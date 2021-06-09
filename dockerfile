FROM rocker/r-ver:4.0.4

ARG WHEN

# required
MAINTAINER Erik Skie <skiex003@umn.edu>

# copy the repo contents into the docker image at `/ensemblemerge`
COPY . /ensemblemerge

# install the dependencies of the R package located at `/portalDS`
RUN apt-get -y update -qq \ 
  && apt-get install -y --no-install-recommends \
    libgsl0-dev \
  && R -e "install.packages('devtools')" \
  && R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')" \
  && R -e "BiocManager::install(c('SummarizedExperiment', 'SingleCellExperiment', 'LoomExperiment'))" \
  && R -e "remotes::install_github('immunogenomics/harmony')" \
  && R -e "remotes::install_version('Seurat', version = '4.0.1')" \
  && R -e "remotes::install_github('satijalab/seurat-wrappers')" \
  && R -e "remotes::install_github('cellgeni/sceasy')" \
  && R -e "remotes::install_github('erikjskie/ensemblemerge', auth_token = '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', force = FALSE, ref = 'development')"
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster.r
\name{Cluster,Seurat,SpectrumCluster-method}
\alias{Cluster,Seurat,SpectrumCluster-method}
\title{Cluster a Seurat object by CIDR (https://github.com/VCCRI/CIDR)}
\usage{
\S4method{Cluster}{Seurat,SpectrumCluster}(x, params, ...)
}
\arguments{
\item{x}{a Seurat object}

\item{params}{a CIDRCluster object}

\item{...}{Additional arguments}
}
\value{
returns a data object with clustering results in meta data
}
\description{
Cluster a Seurat object by CIDR (https://github.com/VCCRI/CIDR)
}
\references{
Lin, P., Troup, M. & Ho, J.W. CIDR: Ultrafast and accurate clustering through imputation for single-cell RNA-seq data. Genome Biol 18, 59 (2017). https://doi.org/10.1186/s13059-017-1188-0
}

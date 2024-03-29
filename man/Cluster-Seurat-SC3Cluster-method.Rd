% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster.r
\name{Cluster,Seurat,SC3Cluster-method}
\alias{Cluster,Seurat,SC3Cluster-method}
\title{Cluster a Seurat object by SC3 (http://bioconductor.org/packages/release/bioc/html/SC3.html)}
\usage{
\S4method{Cluster}{Seurat,SC3Cluster}(x, params, ...)
}
\arguments{
\item{x}{a Seurat object}

\item{params}{a SC3Cluster object}

\item{...}{Additional arguments}
}
\value{
returns a data object with clustering results in meta data
}
\description{
Cluster a Seurat object by SC3 (http://bioconductor.org/packages/release/bioc/html/SC3.html)
}
\references{
Kiselev VY, Kirschner K, Schaub MT, Andrews T, Yiu A, Chandra T, Natarajan KN, Reik W, Barahona M, Green AR, Hemberg M (2017). “SC3 - consensus clustering of single-cell RNA-Seq data.” Nature Methods. http://dx.doi.org/10.1038/nmeth.4236.
}

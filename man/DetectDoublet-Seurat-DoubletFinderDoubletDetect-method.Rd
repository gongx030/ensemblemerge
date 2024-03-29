% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doublet_detect.r
\name{DetectDoublet,Seurat,DoubletFinderDoubletDetect-method}
\alias{DetectDoublet,Seurat,DoubletFinderDoubletDetect-method}
\title{Doublet detection by DoubletFinder (https://github.com/chris-mcginnis-ucsf/DoubletFinder)}
\usage{
\S4method{DetectDoublet}{Seurat,DoubletFinderDoubletDetect}(x, params, ...)
}
\arguments{
\item{x}{a Seurat object}

\item{params}{a DoubletFinderDoubletDetect}

\item{...}{Additional arguments}
}
\value{
returns a data object with doublet removed
}
\description{
Doublet detection by DoubletFinder (https://github.com/chris-mcginnis-ucsf/DoubletFinder)
}
\references{
McGinnis CS, Murrow LM, Gartner ZJ. DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors. Cell Syst. 2019 Apr 24;8(4):329-337.e4. doi: 10.1016/j.cels.2019.03.003. Epub 2019 Apr 3. PMID: 30954475; PMCID: PMC6853612.
}

\name{metaprep-class}
\docType{class}
\alias{metaprep-class}
\title{Class metaprep}
\description{This is a class representation for the effect size estimates and other summary information for a single gene, from a collection of gene expression studies, usually constructed in preparation for a meta-analysis. }
\section{Objects from the Class}{
  Objects can be created using the function \code{metahdep.format} and \code{new}.
}
\section{Slots}{
  \describe{
    \item{\code{theta}:}{Object of class \code{vector} representing the gene's effect size estimates (differential expression measures) from the multiple studies. }
    \item{\code{V}:}{Object of class \code{matrix} representing the sampling variance/covariance matrix of the gene's effect size estimates from the multiple studies. }
    \item{\code{X}:}{Object of class \code{matrix} representing the covariate (or design) matrix for the gene. Covariate information from the multiple studies is represented here. }
    \item{\code{M}:}{Object of class \code{matrix} representing the block diagonal hierarchical structure of effect size estimates from the multiple studies. } 
    \item{\code{max.k}:}{Object of class \code{integer} representing the size of the largest block on the diagonal of \code{M}, i.e., the size of the largest hierarchical dependence group for the gene. }
    \item{\code{row.indices}:}{Object of class \code{matrix} with columns named Study and Row.  Optionally returned by the function \code{metahdep.format()}, to see which gene (Row) in which \code{ES.obj} object (Study) produced the data recorded in each \code{metaprep} object. }
    \item{\code{gene}:}{Object of class \code{character} representing the gene name. }
  }
}
\section{Methods}{
  \describe{
    \item{@}{replace the slot entries }
  }
}

\references{

    Stevens J.R. and Nicholas G. (2009), metahdep: Meta-analysis of hierarchically dependent gene expression studies,
     \emph{Bioinformatics}, 25(19):2619-2620.

    See also the \emph{metahdep} package vignette.

}



\examples{

###
###  See the metahdep package vignette for a full example 
###
data(HGU.prep.list)
HGU.prep.list[[7]]

}

\keyword{classes}

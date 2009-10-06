\name{ES.obj-class}
\docType{class}
\alias{ES.obj-class}
\title{Class ES.obj}
\description{This is a class representation for the effect size estimates and other summary information from a single gene expression study, usually constructed in preparation for meta-analysis. }
\section{Objects from the Class}{
  Objects can be created using the functions \code{getPLM.es} or \code{new}.
}
\section{Slots}{
  \describe{
    \item{\code{gn}:}{Object of class \code{character} representing the probeset IDs of the genes in the study. }
    \item{\code{ES.mat}:}{Object of class \code{matrix} representing the effect size estimates for each gene in the study.
        Rows correspond to probesets and columns correspond to different comparisons or tests of differential expression.  
        If a test of differential expression was performed for different covariate levels, then there will be more than
        one column, so that each row in this matrix represents a vector of effect size estimates for the corresponding probeset
        in the \code{gn} slot. }
    \item{\code{Cov.mat}:}{Object of class \code{matrix}, with each row representing the upper triangle of the variance / covariance
        matrix of the vector of effect size estimates (row in the \code{ES.mat} slot) for the corresponding probeset in the \code{gn} slot.
        Within each row, the order is the same as the result of a call to the \code{upperTriangle(matrix,diag=T)} function in the \emph{gdata} package. }
    \item{\code{chip}:}{Object of class \code{character} representing the chip or array version used in the study. } 
    \item{\code{covariates}:}{Object of class \code{data.frame} representing covariate differences among the columns of the matrix in the \code{ES.mat} slot.
         This object has a row for each column of the matrix in the \code{ES.mat} slot, and a column for each covariate to be considered in the
         meta-analysis, regardless of whether the covariate takes on multiple values in the study represented in this \code{ES.obj} object.
         For best interpretability, columns of the \code{data.frame} in this \code{covariates} slot should be coded as 0/1. }
    \item{\code{dep.grp}}{Object of class \code{integer} representing the dependence group number assigned to the study.
        Studies from the same research team may be considered hierarchically dependendent and share the same \code{dep.grp} value. }
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
data(HGU.DifExp.list)
ES <- HGU.DifExp.list[[1]]
slotNames(ES)
head(ES@gn)
head(ES@ES.mat)
head(ES@Cov.mat)
ES@chip
ES@covariates
ES@dep.grp


}

\keyword{classes}

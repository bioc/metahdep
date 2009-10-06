\name{getPLM.es}
\alias{getPLM.es}
\title{ getPLM.es }
\description{
Calculates effect size estimates for a single study, based on a probe-level model, in preparation for a meta-analysis.  It returns an \code{ES.obj} object containing the result.
}
\usage{
getPLM.es(abatch, trt1, trt2, covariates=NULL, dep.grp=NULL, 
          sub.gn=NULL, bg.norm=TRUE)
}
\arguments{
  \item{abatch}{ An \code{AffyBatch} object containing the data of interest. }
  \item{trt1}{ A vector (or list of vectors) of array indices for treatment level 1 (control).  If more than one test of differential expression is to be performed (for multiple covariate levels, for example), this should be a list of vectors; each \code{trt1} / \code{trt2} vector pair defines a comparison of interest.  }
  \item{trt2}{ A vector (or list of vectors) of array indices for treatment level 2 (treatment).  If more than one test of differential expression is to be performed (for multiple covariate levels, for example), this should be a list of vectors; each \code{trt1} / \code{trt2} vector pair defines a comparison of interest.    }
  \item{covariates}{ (optional) A \code{data.frame} object representing covariate differences, if any, among the comparisons defined by \code{trt1} / \code{trt2} vector pairs. This \code{data.frame} should have a named column for each covariate to be considered in the meta-analysis, regardless of whether the covariate takes on multiple values in the study represented by the \code{abatch} argument.  This \code{data.frame} must have a row for each comparison of interest, as defined by the \code{trt1} / \code{trt2} vector pairs.  Elements of this \code{data.frame} should be coded numerically. }
  \item{dep.grp}{ (optional) A single numeric value representing the dependence group number assigned to the study. Studies from the same research team may be considered hierarchically dependendent and share the same value. }
  \item{sub.gn}{ (optional) A vector of geneNames (probe set ID's); the probe-level model will only be fit for these probesets.  If NULL (default), all probesets are used. }
  \item{bg.norm}{ (optional) A logical value specifying whether or not to perform background correction and normalization before fitting the probe-level model. }
}
\details{
For some subset of probesets in a gene expression study, this function calculates the effect size estimates based on Bolstad's probe-level model (Bolstad 2004), as described in Hu et al. (2006).
Only two-group comparisons (treatment vs. control, for example) are supported.
This is done in preparation for a meta-analysis of multiple gene expression studies.
}
\value{
   An object of class \code{ES.obj}
 }
\author{ John R. Stevens, Gabriel Nicholas }

\references{

    Bolstad B. M. (2004),  \emph{Low-level Analysis of High-density Oligonucleotide Array Data: Background, Normalization and Summarization}, PhD dissertation, U.C. Berkeley.

    Hu P., Greenwood C.M.T., and Beyene J. (2006), Integrative Analysis of Gene Expression Data Including an Assessment of Pathway Enrichment for Predicting Prostate Cancer, 
      \emph{Cancer Informatics} 2006:2 289-300.

    Stevens J.R. and Nicholas G. (2009), metahdep: Meta-analysis of hierarchically dependent gene expression studies,
     \emph{Bioinformatics}, 25(19):2619-2620.

    See also the \emph{metahdep} package vignette.

}


\examples{

###
###  See the metahdep package vignette for a full example 
###


}
\keyword{ models }
\keyword{ htest }
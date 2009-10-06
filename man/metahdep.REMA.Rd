\name{metahdep.REMA}
\alias{metahdep.REMA}
\title{ metahdep.REMA }
\description{
Performs a random effects linear model meta-analysis, allowing for hierarchical dependence.  It returns a list containing the results.
}
\usage{
metahdep.REMA(theta, V, X, M = NULL, dep.groups = NULL, 
              meta.name = "meta-analysis", delta.split = FALSE, 
              center.X = FALSE)
}
\arguments{
  \item{theta}{ A vector of effect size estimates from multiple studies. }
  \item{V}{ The variance/covariance matrix for \code{theta}. Typically, this will be block diagonal (to represent any sampling dependence). }
  \item{X}{ A matrix of covariates for \code{theta}.  At the very least, this must consist of an intercept term.  Other covariates can be included, but there must be more rows than columns in this covariate matrix. }
  \item{M}{ (optional) Used when \code{delta.split=TRUE}.  A block-diagonal matrix describing the hierarchical dependence for the studies (\code{theta}).  One of two ways to specify this is by using the \code{metahdep.format()} function; the other is to use the \code{get.M()} function. }
  \item{dep.groups}{ (optional) Used when \code{delta.split=TRUE}.  A list of vectors/scalars describing the hierarchical dependence groups for the studies (\code{theta}).  This is an alternative to passing an \code{M} matrix. }
  \item{meta.name}{ (optional) A name field for bookkeeping.  This can be any character string.  }
  \item{delta.split}{ (optional) A logical value specifying whether or not to account for hierarchical dependence (i.e., perform delta-splitting).  If \code{TRUE}, then the user needs to pass either a dependence matrix \code{M}, or a \code{dep.groups} list; i.e., one of \code{M} or \code{dep.groups} is REQUIRED when \code{delta.split=TRUE}. }
  \item{center.X}{ (optional) A logical value specifying whether or not to center the columns of \code{X}.  If \code{TRUE}, then the mean from each column will be subtracted from every element in that column (but not for the intercept).  This changes the interpretation of the intercept coefficient estimate from the model fit. }
}
\details{
Takes a vector of effect size estimates, a variance/covariance matrix, and a covariate matrix, and fits a random effects linear model meta-analysis, allowing for hierarchical dependence.
If \code{delta.split=TRUE}, then it performs delta-splitting to account for hierarchical dependence among studies.
When a meta-analysis is to be performed for gene expression data (on a per-gene basis), the \code{metahdep()} function calls this \code{metahdep.REMA()} function for each gene separately.
}
\value{
   A list, with the following named components:
  \item{beta.hats }{A vector of model estimates for the covariates given by \code{X} (it may be a scalar, i.e., a vector of length 1) }
  \item{cov.matrix }{The variance/covariance matrix for the \code{beta.hats} vector }
  \item{beta.hat.p.values }{The [two-sided] p-value(s) for the \code{beta.hats} estimate(s) }
  \item{tau2.hat }{The estimated between-study hierarchical variance tau-square, using the method of moments approach of DerSimonian and Laird. }
  \item{varsigma.hat }{(Only estimated when \code{delta.split=TRUE}.) The estimated within-group hierarchical covariance. }
  \item{Q }{The statistic used to test for model homogeneity / model mis-specification }
  \item{Q.p.value }{The p-value for \code{Q} }
  \item{name }{An optional name field }
}
\author{ John R. Stevens, Gabriel Nicholas }

\references{

  DerSimonian R. and Laird N. (1986), Meta-analysis in clinical trials, \emph{Controlled Clinical Trials}, 7: 177-188.

  Hedges L. V. and Olkin I (1985), \emph{Statistical methods for meta-analysis}, San Diego, CA: Academic Press.

  Stevens J.R. and Doerge R.W. (2005), Combining Affymetrix Microarray Results, \emph{BMC Bioinformatics}, 6:57.

  Stevens J.R. and Taylor A.M. (2009), Hierarchical Dependence in Meta-Analysis, 
   \emph{Journal of Educational and Behavioral Statistics}, 34(1):46-73.

  See also the \emph{metahdep} package vignette.

}


\examples{

###
### Example 1: gene expression data
### - this uses one gene from the HGU.prep.list object

# load data and extract components for meta-analysis (for one gene)
data(HGU.prep.list)
gene.data <- HGU.prep.list[[7]]
theta <- gene.data@theta
V <- gene.data@V
X <- gene.data@X
M <- gene.data@M
dep.grps <- list(c(1:2),c(4:6))
gene.name <- gene.data@gene

# fit a regular REMA (no hierarchical dependence)
results <- metahdep.REMA(theta, V, X, meta.name=gene.name)
results

# fit hierarchical dependence model (with delta-splitting), 
# using two different methods for specifying the dependence structure
results.dsplitM <- metahdep.REMA(theta, V, X, M, delta.split=TRUE, 
  meta.name=gene.name, center.X=TRUE)
results.dsplitM
results.dsplitd <- metahdep.REMA(theta, V, X, dep.groups=dep.grps, 
  delta.split=TRUE, meta.name=gene.name, center.X=TRUE)
results.dsplitd


###
###  Example 2: glossing data
###  - this produces part of Table 6 in the Stevens and Taylor JEBS paper.

data(gloss)
dep.groups <- list(c(2,3,4,5),c(10,11,12))
REMA.ds <- metahdep.REMA(gloss.theta, gloss.V, gloss.X, center.X=TRUE, 
  delta.split=TRUE, dep.groups=dep.groups)
round(cbind(t(REMA.ds$beta.hats), sqrt(diag(REMA.ds$cov.matrix)), 
  t(REMA.ds$beta.hat.p.values)),4)



}
\keyword{ models }
\keyword{ htest }
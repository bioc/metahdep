\name{metahdep.HBLM}
\alias{metahdep.HBLM}
\title{metahdep.HBLM }
\description{
Performs a meta-analysis by fitting a hierarchical Bayes linear model, allowing for hierarchical dependence.
}
\usage{
metahdep.HBLM(theta, V, X, M = NULL, dep.groups = NULL, 
              meta.name = "meta-analysis", center.X = FALSE, 
              delta.split = FALSE, n = 10, m = 10, 
              two.sided = FALSE)
}
\arguments{
  \item{theta}{ A vector of effect size estimates from multiple studies. }
  \item{V}{ The variance/covariance matrix for \code{theta}. Typically, this will be block diagonal (to represent any sampling dependence). }
  \item{X}{ A matrix of covariates for \code{theta}.  At the very least, this must consist of an intercept term.  Other covariates can be included, but there must be more rows than columns in this covariate matrix. }
  \item{M}{ (optional) Used when \code{delta.split=TRUE}.  A block-diagonal matrix describing the hierarchical dependence for the studies (\code{theta}).  One of two ways to specify this is by using the \code{metahdep.format()} function; the other is to use the \code{get.M()} function. }
  \item{dep.groups}{ (optional) Used when \code{delta.split=TRUE}.  A list of vectors/scalars describing the hierarchical dependence groups for the studies (\code{theta}).  This is an alternative to passing an \code{M} matrix. }
  \item{meta.name}{ (optional) A name field for bookkeeping.  This can be any character string.  }
  \item{center.X}{ (optional) A logical value specifying whether or not to center the columns of \code{X}.  If \code{TRUE}, then the mean from each column will be subtracted from every element in that column (but not for the intercept).  This changes the interpretation of the intercept coefficient estimate from the model fit. }
  \item{delta.split}{ (optional) A logical value specifying whether or not to account for hierarchical dependence (i.e., perform delta-splitting).  If \code{TRUE}, then the user needs to pass either a dependence matrix \code{M}, or a \code{dep.groups} list; i.e., one of \code{M} or \code{dep.groups} is REQUIRED when \code{delta.split=TRUE}. }
  \item{n}{ (optional) An even integer telling how many steps to use when doing the numerical integration over tau, the square root of the between-study hierarchical variance.  The integration is done on the log-logistic prior, split into the 4 quartiles.  This number n specifies how many steps to take within each quartile. }
  \item{m}{ (optional) An even integer telling how many steps to use when doing the numerical integration over varsigma (given tau), the between-study hierarchical covariance.  This is only used when \code{delta.split=TRUE}. The integration is done on the uniform prior, for each value of tau. }
  \item{two.sided}{ (optional) A logical value to determine whether to return the 2-sided p-values or default [one-sided positive] posterior probabilities for the parameter estimates. }
}
\details{
Takes a vector of effect size estimates, a variance/covariance matrix, and a covariate matrix, and fits a hierarchical Bayes linear model.  If \code{delta.split=TRUE}, then it performs delta-splitting to account for hierarchical dependence among studies.  The main parameters (beta) are given normal priors, the square root of the hierarchical variance (tau) is given a log-logistic prior, and the hierarchical covariance (varsigma) is given a uniform prior; see the Stevens and Taylor reference for details.
When a meta-analysis is to be performed for gene expression data (on a per-gene basis), the \code{metahdep()} function calls this \code{metahdep.HBLM} function for each gene separately.
}
\value{
   A list, with the following named components:
  \item{beta.hats }{A vector of model estimates for the covariates given by \code{X} (it may be length 1 i.e. scalar) }
  \item{cov.matrix }{The variance/covariance matrix for the \code{beta.hats} vector }
  \item{beta.hat.p.values }{The p-value(s) for the \code{beta.hats} estimate(s) }
  \item{tau.hat }{The posterior mean for tau (not tau-square).  An estimate for tau-square is E(square(tau) [given data]) = tau.var + square(tau.hat)  }
  \item{tau.var }{The posterior variance for tau (not tau-square). }
  \item{varsigma.hat }{The posterior mean for varsigma. }
  \item{varsigma.var }{The posterior variance for varsigma. }
  \item{tau.varsigma.cov }{The posterior covariance for tau and varsigma. }
  \item{name }{An optional name field }
}
\author{ John R. Stevens, Gabriel Nicholas }


\references{

DuMouchel W. H. and Harris J. H. (1983), Bayes methods for combining the results of
cancer studies in humans and other species, \emph{Journal of the American Statistical Association},
78(382), 293-308.

DuMouchel W.H. and Normand S.-L. (2000), Computer-modeling and graphical strategies
for meta-analysis, in D. K. Stangl and D. A. Berry (Eds.), \emph{Meta-analysis in medicine
and health policy}, pp. 127-178. New York: Marcel Dekker.

Stevens J.R. and Doerge R.W. (2005), A Bayesian and Covariate Approach to Combine Results
   from Multiple Microarray Studies, \emph{Proceedings of Conference on Applied Statistics in
   Agriculture}, pp. 133-147.

Stevens J.R. and Nicholas G. (2009), metahdep: Meta-analysis of hierarchically dependent gene expression studies,
     \emph{Bioinformatics}, 25(19):2619-2620.

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

# fit a regular HBLM (no hierarchical dependence)
results <- metahdep.HBLM(theta, V, X, meta.name=gene.name, 
  center.X=TRUE, two.sided=TRUE)
results

# fit hierarchical dependence model (with delta-splitting), 
# using two different methods for specifying the dependence structure
results.dsplitM <- metahdep.HBLM(theta, V, X, M, delta.split=TRUE, 
  meta.name=gene.name, center.X=TRUE, two.sided=TRUE)
results.dsplitM
results.dsplitd <- metahdep.HBLM(theta, V, X, dep.groups=dep.grps, 
  delta.split=TRUE, meta.name=gene.name, center.X=TRUE, two.sided=TRUE)
results.dsplitd


###
###  Example 2: glossing data
###  - this produces part of Table 5 in the Stevens and Taylor JEBS paper.

data(gloss)
dep.groups <- list(c(2,3,4,5),c(10,11,12))
HBLM.ds <- metahdep.HBLM(gloss.theta, gloss.V, gloss.X, center.X=TRUE, 
  two.sided=TRUE, delta.split=TRUE, dep.groups=dep.groups, n=20, m=20)
round(cbind(HBLM.ds$beta.hats, HBLM.ds$beta.hat.p.values),4)


}
\keyword{ models }
\keyword{ htest }

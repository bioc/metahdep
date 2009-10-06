\name{metahdep.FEMA}
\alias{metahdep.FEMA}
\title{ metahdep.FEMA }
\description{
Performs a fixed effects linear model meta-analysis.  It returns a list containing the results.
}
\usage{
metahdep.FEMA(theta, V, X, meta.name = "meta-analysis", 
              center.X = FALSE)
}
\arguments{
  \item{theta}{ A vector of effect size estimates from multiple studies. }
  \item{V}{ The variance/covariance matrix for \code{theta}. Typically, this will be block diagonal (to represent any sampling dependence). }
  \item{X}{ A matrix of covariates for \code{theta}.  At the very least, this must consist of an intercept term.  Other covariates can be included, but there must be more rows than columns in this covariate matrix. }
  \item{meta.name}{ (optional) A name field for bookkeeping.  This can be any character string. }
  \item{center.X}{ (optional) A logical value specifying whether or not to center the columns of \code{X}.  If \code{TRUE}, then the mean from each column will be subtracted from every element in that column (but not for the intercept).  This changes the interpretation of the intercept coefficient estimate from the model fit. }
}
\details{
Takes a vector of effect size estimates, a variance/covariance matrix, and a covariate matrix, and fits a fixed effects linear model meta-analysis.
When a meta-analysis is to be performed for gene expression data (on a per-gene basis), the \code{metahdep()} function calls this function for each gene separately.
}
\value{ A list with the following named components:
  \item{beta.hats }{A vector of model estimates for the covariates given by \code{X} (it may be a scalar, i.e., a vector of length 1 )}
  \item{cov.matrix }{The variance/covariance matrix for the \code{beta.hats} estimate(s)}
  \item{beta.hat.p.values }{The [two-sided] p-value(s) for the \code{beta.hats} estimate(s)}
  \item{Q }{The statistic used to test for model homogeneity / model mis-specification}
  \item{Q.p.value }{The p-value for \code{Q}}
  \item{name }{An optional name field}
}
\author{ John R. Stevens, Gabriel Nicholas }

\references{

  Hedges L. V. and Olkin I (1985), \emph{Statistical methods for meta-analysis}, San Diego, CA: Academic Press.

  Stevens J.R. and Doerge R.W. (2005), Combining Affymetrix Microarray Results, \emph{BMC Bioinformatics} 6:57.

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
gene.name <- gene.data@gene

# fit a regular FEMA (no hierarchical dependence)
results <- metahdep.FEMA(theta, V, X, meta.name=gene.name, center.X=TRUE)
results

###
###  Example 2: glossing data
###  - this produces part of Table 5 in the Stevens and Taylor JEBS paper.

data(gloss)
FEMA <- metahdep.FEMA(gloss.theta, gloss.V, gloss.X, center.X=TRUE)
round(cbind( t(FEMA$beta.hats), t(FEMA$beta.hat.p.values)),4)


}
\keyword{ models }
\keyword{ htest }

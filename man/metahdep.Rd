\name{metahdep}
\alias{metahdep}
\title{ metahdep }
\description{
Takes a list of \code{metaprep} objects and performs the specified meta-analysis on each element.  Intended mainly for meta-analyzing the results of multiple gene expression studies. }
\usage{
metahdep(prep.list, genelist = NULL, method = "HBLM", n = 10, m = 10, 
         center.X = FALSE, delta.split = FALSE, return.list = FALSE, 
         two.sided = TRUE)
}
\arguments{
  \item{prep.list}{ A list of \code{metaprep} class objects as returned by the \code{metahdep.format()} function. }
  \item{genelist}{ (optional) A subsetting parameter.  A vector of gene/probeset names on which to perform the meta-analyses. }
  \item{method}{ (optional) One of: "FEMA" - fixed effects meta-analysis, "REMA" - random effects meta-analysis, or "HBLM" - hierarchical Bayes linear model.  This defaults to "HBLM". }
  \item{n}{ (optional) An even integer specifying the number of steps to take over each quartile in the numerical integration over tau when doing HBLM.  See \code{metahdep.HBLM}. }
  \item{m}{ (optional) An even integer specifying the number of steps to take in the numerical integration over varsigma (given tau) when doing HBLM.  See \code{metahdep.HBLM}. }
  \item{center.X}{ (optional) A logical value specifying whether or not to center the columns of the covariate matrices.  If \code{TRUE}, then for the covariate matrix of each \code{metaprep} object, the mean each non-intercept column will be subtracted from every element in the column prior to the meta-analysis.  This changes the interpretation of the intercept coefficient estimate from the model fit. }
  \item{delta.split}{ (optional) A logical value specifying whether or not to account for hierarchical dependence via delta-splitting.  Only used in methods "REMA" and "HBLM".  If \code{TRUE}, then each \code{metaprep} object must include a dependence matrix M. }
  \item{return.list}{ (optional) A logical value specifying whether to return the results as a list of lists rather than as a \code{data.frame}. The default is \code{FALSE}.}
  \item{two.sided}{ (optional) A logical value specifying whether to transform the posterior probabilities from the HBLM method.  The default \code{TRUE} returns 2-sided p-values for the parameter estimates for convenience in interpretation.  If this is set to FALSE, then it will return 1-sided posterior probabilities representing P( beta[j] > 0 | data ). }
}
\value{
  Returns a \code{data.frame} by default.  The exact contents of the \code{data.frame} will vary depending on the \code{method} argument.  The row names of the \code{data.frame} will be the gene names from the \code{prep.list} argument.  For all \code{method} options, the first several columns of the resulting \code{data.frame} will be the model parameter estimates (beta hats).  The next group of columns will be the elements of the variance/covariance matrix for the beta hats.  The next group of columns will be the p-values for the parameter estimates.  The remaining columns will change depending on the method. 

  For FEMA (and REMA), the remaining columns are the Q statistic and its p-value -- testing for model homogeneity.

  For HBLM, the remaining columns are the posterior mean and variance of tau, the posterior mean and variance of varsigma, and the posterior covariance of tau and varsigma.

  All columns in the \code{data.frame} have meaningful names to aid their interpretation.

}


\author{ John R. Stevens, Gabriel Nicholas }


\references{

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
data(HGU.prep.list)

##  do FEMA and REMA, and view the results
FEMA.results <- metahdep(HGU.prep.list, method="FEMA", center.X=TRUE)
head(FEMA.results)

REMA.results <- metahdep(HGU.prep.list, method="REMA", center.X=TRUE)
head(REMA.results)

##  get a small subset of genes
##  some of these may not be suitable for all methods 
##  (there may not be enough data for that gene)
data(HGU.newnames)
set.seed(123)
gene.subset <- sample(HGU.newnames$new.name, 50)

##  view results from REMA and HBLM with delta splitting on subset of genes
REMA.dsplit.results <- metahdep(HGU.prep.list, method="REMA", 
  genelist=gene.subset, delta.split=TRUE, center.X=TRUE)
head(REMA.dsplit.results)

HBLM.dsplit.results <- metahdep(HGU.prep.list, method="HBLM", 
  genelist=gene.subset, delta.split=TRUE, center.X=TRUE)
head(HBLM.dsplit.results)
}
\keyword{ models }
\keyword{ htest }


\name{gloss}
\docType{data}
\alias{gloss}
\alias{gloss.Table1}
\alias{gloss.X}
\alias{gloss.theta}
\alias{gloss.V}
\title{gloss:  Data from the JEBS glossing paper }
\description{
  This includes the following four objects:
\tabular{ll}{
  \tab \cr
  \code{gloss.Table1} \tab a data.frame containing the contents of Table 1 in the paper \cr
  \tab (sample sizes, sample means, sample SDs, and covariate information from each study) \cr
  \code{gloss.X} \tab a matrix representing the design matrix X for the meta-analysis in the paper \cr
  \code{gloss.theta} \tab a vector, representing the effect size estimates as summarized in Table 3 of the paper \cr
  \code{gloss.V} \tab a matrix, representing the covariance matrix of effect size estimates, including \cr
  \tab sampling dependence in off-diagonal elements, as summarized in Table 3 of paper \cr
  }
}
\usage{data(gloss)}
\format{
  This object contains the four objects described above.
}
\details{

This data set summarizes 13 experiments with 18 study reports, all involving the effect of native-language (L1) vocabulary aids on second language (L2) reading comprehension. Some experiments produced multiple study reports, creating a dependence structure among the resulting effect size estimates.

The intended use for these data is to demonstrate the methods coded in the \emph{metahdep} package.
}
\references{ Stevens J.R. and Taylor A.M. Hierarchical Dependence in Meta-Analysis. To appear, Journal of Educational and Behavioral Statistics. }
\examples{
data(gloss)
# Look at Table 1
gloss.Table1
}
\keyword{datasets}

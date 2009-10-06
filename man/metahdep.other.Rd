\name{metahdep.other}
\alias{metahdep.list2dataframe}
\alias{LinMod.MetAn.dep.REMA}
\alias{LinMod.REMA.dep}
\alias{LinMod.REMA.delta.split}
\alias{LinMod.HBLM.fast.dep}
\alias{new.LinMod.HBLM.fast.dep.delta.split}
\alias{LinMod.MetAn.dep.FEMA}
\alias{metahdep.check.X}
\alias{get.M}
\alias{tr}
\alias{id}
\alias{center.columns}
\alias{mod}
\alias{get.varsigma.v}
\title{ metahdep.other }
\description{
Miscellaneous functions used internally by the \emph{metahdep} package's main functions (\code{metahdep}, \code{metahdep.FEMA}, \code{metahdep.REMA}, \code{metahdep.HBLM}, and \code{metahdep.format}):
\tabular{ll}{
\tab \cr
\code{metahdep.list2dataframe}  \tab convert list to \code{data.frame} \cr
\code{LinMod.MetAn.dep.REMA}  \tab  REMA meta-analysis  \cr
\code{LinMod.REMA.dep}  \tab used by \code{LinMod.MetAn.dep.REMA} to estimate parameters  \cr
\code{LinMod.REMA.delta.split}    \tab  REMA (with delta-splitting)  \cr
\code{LinMod.HBLM.fast.dep}    \tab HBLM (no delta-splitting)  \cr
\code{new.LinMod.HBLM.fast.dep.delta.split}    \tab  HBLM (with delta-splitting)  \cr
\code{LinMod.MetAn.dep.FEMA}    \tab  FEMA  \cr
\code{metahdep.check.X}   \tab  check design matrix X, and drop columns if necessary \cr
\tab to make full rank  \cr
\code{get.M}  \tab create block diagonal M matrix, given dependence structure  \cr
\code{tr}  \tab calculate trace of matrix   \cr
\code{id}  \tab create identity matrix   \cr
\code{center.columns}  \tab center all non-intercept columns of design matrix X   \cr
\code{mod}  \tab  mod function   \cr
\code{get.varsigma.v}  \tab get varsigma values for HBLM delta-splitting model \cr
\tab \cr
  }
}

\usage{

}
\arguments{

}
\value{

}
\author{ John R. Stevens, Gabriel Nicholas }


\references{

    Stevens J.R. and Nicholas G. (2009), metahdep: Meta-analysis of hierarchically dependent gene expression studies,
     \emph{Bioinformatics}, 25(19):2619-2620.

   Stevens J.R. and Taylor A.M. (2009), Hierarchical Dependence in Meta-Analysis, 
     \emph{Journal of Educational and Behavioral Statistics}, 34(1):46-73.

    See also the \emph{metahdep} package vignette.

}



\examples{
## Create the M matrix for the glossing example
## - here, studies 2-5 are one hierarchically dependent group (Baumann),
## and studies 10-12 are another hierarchically dependent group (Joyce)
data(gloss)
dep.groups <- list(c(2:5),c(10:12))
M <- get.M(length(gloss.theta),dep.groups)

}
\keyword{ models }

\name{HGU.prep.list}
\alias{HGU.prep.list}
\docType{data}
\title{HGU.prep.list}
\description{
An illustrative example of a list, where each element is a \code{metaprep} class object with data for a particular gene.  
It comes from the study data held within the \code{HGU.DifExp.list} object.
}
\usage{data(HGU.prep.list)}
\format{
  This object is a list of \code{metaprep} objects.
}
\details{
Each element of \code{HGU.prep.list} is a gene-specific object of class \code{metaprep}.

\code{HGU.prep.list} is an example of an object created by the \code{metahdep.format()} function; see the help file for the \code{metahdep.format()} function and the \emph{metahdep} package vignette for details on its construction.
\code{HGU.prep.list} is an example object used as an argument to the \code{metahdep()} function.  
The individual elements of this object can be extracted and passed as arguments to the more general meta-analysis functions, \code{metahdep.HBLM()}, \code{metahdep.REMA()}, and \code{metahdep.FEMA()}.
}

\references{

    Stevens J.R. and Nicholas G. (2009), metahdep: Meta-analysis of hierarchically dependent gene expression studies,
     \emph{Bioinformatics}, 25(19):2619-2620.

    See also the \emph{metahdep} package vignette.

}



\examples{
data(HGU.prep.list)

HGU.prep.list[[7]]
## etc...
}
\keyword{datasets}

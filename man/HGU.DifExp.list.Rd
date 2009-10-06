\name{HGU.DifExp.list}
\docType{data}
\alias{HGU.DifExp.list}
\title{HGU.DifExp.list:  A list of 4 ES.obj objects }
\description{
An illustrative example of a list summarizing several studies of gene expression.  This is used as an example in the meta-analyis of hierarchically dependent gene expression studies.
}
\usage{data(HGU.DifExp.list)}
\format{
  This object is a list containing 4 \code{ES.obj} objects.  
  Each \code{ES.obj} object represents the results from a separate gene expression study.

}
\details{
This object has been assembled from existing data as an artificial example; see the vignette for details on its construction.
In this example, four studies were conducted, and can be summarized as follows:
\tabular{cccc}{
      Study \tab   Lab  \tab Tissue \tab  Chip  \cr
        1   \tab    1   \tab    0   \tab hgu133a \cr
        1   \tab    1   \tab    1   \tab hgu133a \cr
        2   \tab    1   \tab    0   \tab hgu95a \cr
        3   \tab    2   \tab    0   \tab hgu95av2 \cr 
        4   \tab    3   \tab    1   \tab hgu133b \cr
  }
Notice that study 1 involved two tissue types.  The vignette shows how this example supposes that sampling dependence was 
introduced in study 1 by fitting a gene-specific model with both tissue types simultaneously.
Hierarchical dependence is also present in these data because studies 1 and 2 were conducted by the same lab.
Each element of \code{HGU.DifExp.list} is an \code{ES.obj} object in the same format as returned by the \code{getPLM.es()} function.
Look at the elements of the list (and the vignette) to get an idea of how the data should be laid out.

The intended use for these data is to demonstrate a meta-analysis procedure that accounts for hierarchical dependence between studies.  
The idea is that results from different studies from the same lab might be dependent.  
This is an example object that is to be passed as an argument to the \code{metahdep.format()} function.
}

\references{

    Stevens J.R. and Nicholas G. (2009), metahdep: Meta-analysis of hierarchically dependent gene expression studies,
     \emph{Bioinformatics}, 25(19):2619-2620.

    See also the \emph{metahdep} package vignette.

}


\examples{
data(HGU.DifExp.list)
head(HGU.DifExp.list[[1]]@ES.mat)
HGU.DifExp.list[[1]]@covariates
##  etc.
}
\keyword{datasets}

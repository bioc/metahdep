\name{HGU.newnames}
\alias{HGU.newnames}
\docType{data}
\title{HGU.newnames}
\description{
An illustrative example of a \code{data.frame} with 3 columns.  
The first column is named "chip", the second "old.name", and the third "new.name".  
The rows each hold the name of a chip type, a chip-specific probeset name, and a common name used to match probesets across different chip versions.
}
\usage{data(HGU.newnames)}
\format{
  A \code{data.frame} with observations on the following 3 variables for a subset of probesets on different chip types.
  \describe{
    \item{\code{chip}}{a character vector specifying the chip type}
    \item{\code{old.name}}{a character vector specifying the probeset name on the chip type}
    \item{\code{new.name}}{a character vector specifying the common identifier, such as an Entrez Gene ID, for the probeset on the chip type}
  }
}
\details{
This is an example of a \code{newnames} argument that is required by the \code{metahdep.format} function.  
When paired with a list of \code{ES.obj} class objects (see \code{HGU.DifExp.list}) this allows the \code{metahdep.format()} 
function to assemble all of the information from all of the studies for a specific gene.  
The \code{new.name} is a 'common' identifier, e.g., an Entrez Gene ID.  Different studies may use different chip types, or different versions of chips, where information for a gene with a particular Entrez Gene ID may have a different probeset name on each chip type.  
This \code{newnames} argument is meant to facilitate the matching of gene information across different chip types.
See the \emph{metahdep} package vignette for more details on the construction of this object.
}
\examples{
data(HGU.newnames)
head(HGU.newnames)
}
\keyword{datasets}

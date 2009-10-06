\name{metahdep.format}
\alias{metahdep.format}
\title{ metahdep.format }
\description{
This function is intended to facilitate the meta-analysis of multiple gene expression studies.  
This function takes the results from a number of studies and collects together the data for each gene, in preparation for a meta-analysis. }
\usage{
metahdep.format(ES.obj.list, newnames, min.var = 0.0001, 
                include.row.indices = FALSE, show.warnings = FALSE, 
                pd.verify = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ES.obj.list}{ A list object containing the results from multiple gene expression studies.  
                      Each element of the list is an object of class \code{ES.obj}. }
  \item{newnames}{ A \code{data.frame} object that describes how to merge the results from different studies.  It must have 3 columns, and a row for each probeset to be considered.  The first column must be named "chip" and contain a name of a chip version used in one of the studies, corresponding to the "chip" slots of the \code{ES.obj.list} argument.  The second column must be named "old.name" and must hold the chip-specific name for each probeset.  The third column must be named "new.name" and should hold the general name (Entrez ID, for example) that can be used to identify common probesets on the different chip versions.  See the \code{HGU.newnames} object in the given example code (and its creation in the vignette) for guidance on how to build this object. }
  \item{min.var}{ (optional) A positive real number that acts as a lower bound on the allowed variances for any measure of differential expression.  This might be used to guard against specious claims of significance due to naturally low variance. }
  \item{include.row.indices}{ (optional) A logical value to determine whether or not to include the study and row indices for the data for each gene in the returned list of \code{metaprep} objects.  Only used for debugging purposes, and so defaults to FALSE. }
  \item{show.warnings}{ (optional) A logical value to determine whether or not to display warnings in certain situations, like if a gene is expected in a particular study but not found, or if a gene is not found in any study, or on some other instances.  This may sometimes cause a large number of uninteresting warnings to be displayed, and so defaults to FALSE. }
  \item{pd.verify}{ (optional) A logical value to determine whether or not to check the generated variance/covariance matrix for positive-definiteness.  Since this should always be the case, this would indicate a problem somewhere in the data -- more specifically, in the covariance values of one of the studies in the \code{ES.obj.list} argument.  Used primarily for debugging, and so defaults to FALSE. }
}
\details{
Each element of the returned list is a \code{metaprep} object summarizing effect size data for a single gene.  This list of \code{metaprep} gene information is passed to the \code{metahdep()} function for meta-analysis. 
}
\value{
  A list, where each element is a gene-specific object of class \code{metaprep}.
}
\author{ John R. Stevens, Gabriel Nicholas }

\references{

    Stevens J.R. and Nicholas G. (2009), metahdep: Meta-analysis of hierarchically dependent gene expression studies,
     \emph{Bioinformatics}, 25(19):2619-2620.

    See also the \emph{metahdep} package vignette.

}



\examples{
##  load a pre-made list of ES.obj's and newnames data.frame.  These objects hold
##  data in a format suitable for use in the metahdep.format function.
data(HGU.DifExp.list)
data(HGU.newnames)

##  now call the format function;
##  this may take anywhere from several seconds to several minutes, 
##  depending on the speed of the computer and the number of genes under
##  consideration
HGU.prep.list <- metahdep.format(HGU.DifExp.list, HGU.newnames)
}
\keyword{ models }

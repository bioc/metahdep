# Copyright 2008 John R. Stevens
# Distributed as part of the metahdep package, under the terms of the GNU General Public License (see DESCRIPTION file)

`metahdep` <-
function(prep.list, genelist=NULL, method="HBLM", n=10, m=10, center.X=FALSE, delta.split=FALSE, return.list=FALSE, two.sided=TRUE)
{
  ######################################################################
  ##  INITIALIZE
  ######################################################################

  ##  check arguments
  if (!is.element(method, c("FEMA", "REMA", "HBLM")))
  {
    methods.string <- "\"FEMA\" \"REMA\" \"HBLM\""
    cat("Invalid method: ", method, "\n")
    cat("Valid methods are: ", methods.string, "\n")
    return(NULL)
  }
  if (method=="HBLM")
  {
    if(((mod(n,2)==1) | n < 2) | (delta.split & ((mod(m,2)==1) | m < 2)))
    {
      cat("Error -- Numerical integration parameters n and m must both be even and >= 2\n")
      return(NULL)
    }
    else if((mod(n,2)==1) | n < 2)
    {
      cat("Error -- Numerical integration parameter n must be even and >= 2\n")
      return(NULL)
    }
  }

  ##  need to get the vector of indices into the prep list corresponding to the elements of the passed genelist subset
  ##  also, find the maximum rank of the design matrices to be able to make the results data frame the correct size
  ##  if the genelist is NULL then no genelist was passed, and so the full prep list should be analyzed
  if (!is.null(genelist))
    genelist <- unique(genelist)
  metah.update.indices <- trunc(seq(1, length(prep.list), length.out=22))
  cat("|----Initializing----|\n")
  max.rank <- 1
  index.vector <- rep(NA, length(prep.list))
  cov.names <- names(prep.list[[1]]@X)
  for (i in 1:length(prep.list))
  {
    temp.element <- prep.list[[i]]
    if (is.element(i, metah.update.indices))
    {
      cat("*")
      flush.console()
    }
    ##  only include the index for the gene if it is present in the prep list AND
    ##  has more than a single differential expression measure
    if (!is.null(genelist))
    {
      if (is.element(temp.element@gene, genelist) & length(temp.element@theta) > 1)
        index.vector[i] <- i
    }
    else if (length(temp.element@theta) > 1)
      index.vector[i] <- i

    ##  check the covariate matrix for square-ness
    ##  if it's square, then cut off the rightmost column
    ##  this could be moved into the prep.list format function
    current.X <- as.matrix(temp.element@X)
    if (ncol(current.X) == nrow(current.X))
    {
      current.X <- current.X[,1:(ncol(current.X)-1)]
      prep.list[[i]]@X <- as.matrix(current.X)
    }
    current.rank <- qr(current.X)$rank
    if (current.rank > max.rank)
    {
      max.rank <- current.rank
      cov.names <- colnames(current.X)
    }
  }
  cat("\n")

  index.vector <- index.vector[!is.na(index.vector)]
  num.genes <- length(index.vector)
  if (num.genes < 1)
  {
    cat("Warning:  No genes were found with enough data to perform a meta-analysis\n")
    return(NULL)
  }

  output.list <- list()
  metah.update.indices <- trunc(seq(1, length(index.vector), length.out=22))
  if (method == "REMA")
  {
    if (delta.split)
      cat("\nFitting Random Effects model (with delta splitting)")
    else
      cat("\nFitting Random Effects model")
  }
  else   if (method == "HBLM")
  {
    if (delta.split)
      cat("\nFitting Hierarchical Bayes Linear Model (with delta splitting)")
    else
      cat("\nFitting Hierarchical Bayes Linear Model")
  }
  else
    cat("\nFitting Fixed Effects model")
  cat("\n|-------Working------|\n")

  ######################################################################
  ##  MAIN LOOP
  ######################################################################

  for (i in 1:num.genes)
  {
    if (is.element(i, metah.update.indices))
    {
      cat("*")
      flush.console()
    }

    current.gene <- prep.list[[index.vector[i]]]

    if(center.X & ncol(current.gene@X) > 1)
      current.gene@X <- center.columns(current.gene@X)

    max.k <- current.gene@max.k
    if (method == "REMA")
    {
      ##  REMA returns:  beta.hats, cov.matrix, beta.hat.p.values, sigma.hat, varsigma.hat, (Q), (Q.pvalue), gene
      ##  only allow delta splitting if it makes sense (i.e. nonzero M matrix <==> max.k > 0)
      if ((delta.split==TRUE) & (max.k > 0))
        output.list[[i]] <- LinMod.REMA.delta.split(current.gene)
      else
        output.list[[i]] <- LinMod.MetAn.dep.REMA(current.gene)
    }
    else if (method == "HBLM")
    {
      ##  only allow delta splitting if it makes sense (i.e. nonzero M matrix, and also max.k > 0)
      if (!delta.split | max.k == 0)
        output.list[[i]] <- LinMod.HBLM.fast.dep(current.gene, n=n, two.sided=two.sided)
      else
        output.list[[i]] <- new.LinMod.HBLM.fast.dep.delta.split(current.gene, n=n, m=m, two.sided=two.sided)
    }
    else ## method == "FEMA"
    {
      ##  FEMA function returns:  beta.hats, cov.matrix, Q, Q.p.value, beta.hat.p.values, gene
      ##  OR it returns NULL if there was a problem with the gene data (e.g. X not full rank)
      output.list[[i]] <- LinMod.MetAn.dep.FEMA(current.gene)
    }
  }
  cat("\n")
  cat(num.genes, "unique genes were found with enough data to perform a meta-analysis\n")

  ##  if return.list=TRUE then return the results as a list
  if (return.list)
    return(output.list)

  ##  otherwise, put the results into a data frame
  output.frame <- metahdep.list2dataframe(output.list, max.rank=max.rank, method=method, cov.names=cov.names)
  return(output.frame)
}


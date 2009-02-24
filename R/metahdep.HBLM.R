# Copyright 2008 John R. Stevens
# Distributed as part of the metahdep package, under the terms of the GNU General Public License (see DESCRIPTION file)

`metahdep.HBLM` <-
function(theta, V, X, M=NULL, dep.groups=NULL, meta.name="meta-analysis", center.X=FALSE, delta.split=FALSE, n=10, m=10, two.sided=FALSE)
{
  max.k <- 1
  if (center.X)
    X <- center.columns(X)
  X <- metahdep.check.X(X)
  if (is.null(X))
    return(NULL)

  ##  test the arguments
  if (delta.split)
  {
    if (n %% 2 > 0 | m %% 2 > 0 | n < 2 | m < 2)
    {
      cat("Error:  Numerical integration parameters n and m must both be even and greater than 1\n")
      return(NULL)
    }

    if (is.null(M) & is.null(dep.groups))
    {
      cat("Error:  Need a dependence matrix (M) or dependence group information (dep.groups) to perform delta splitting.\n")
      return(NULL)
    }

    if (is.null(M))
      M <- get.M(length(theta), dep.groups)
    if (is.null(M))
      return(NULL)

    if (sum(M) == 0)
    {
      cat("Error:  Unable to perform delta splitting with the current data\n")
      cat("At least 2 observations must come from the same dependence group, and\n")
      cat("there must also be data from at least 2 different dependence groups\n")
      return(NULL)
    }

    if (sum(!is.element(M, c(0,1))) > 0 | !isSymmetric(M) | sum(diag(M)) > 0)
    {
      cat("Error:  Improperly constructed dependence matrix (M)\n")
      cat("M must be symmetric, block diagonal, consist of 0/1, with 0 on the diagonal\n")
      return(NULL)
    }

    ##  find max.k
    max.k <- max(apply(M,1,sum))+1

    if( length(unique(c(length(theta),nrow(V),ncol(V),nrow(X),nrow(M),ncol(M)))) > 1 )
    {
      cat("Error:  Passed arguments are not conformable\n")
      return(NULL)
    }
  }
  if (!delta.split)
  {
    if (n %% 2 > 0 | n < 2)
    {
      cat("Error:  Numerical integration parameter n must be even and greater than 1\n")
      return(NULL)
    }
    if( length(unique(c(length(theta),nrow(V),ncol(V),nrow(X)))) > 1 )
    {
      cat("Error:  Passed arguments are not conformable\n")
      return(NULL)
    }
    test.V <- try(chol(V))
    if (!is.matrix(test.V))
    {
      cat("Error:  The variance/covariance matrix V is not positive definite\n")
      return(NULL)
    }

    if(is.null(M)){M <- 0*V}
  }

  metahdep.info <- new("metaprep", theta=theta, V=V, X=X, M=M, max.k=as.integer(max.k), gene=meta.name)

  if (delta.split)
    return.list <- new.LinMod.HBLM.fast.dep.delta.split(metahdep.info, n=n, m=m, two.sided=two.sided)
  else
    return.list <- LinMod.HBLM.fast.dep(metahdep.info, n=n, two.sided=two.sided)

  return(return.list)
}


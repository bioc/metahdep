# Copyright 2008 John R. Stevens
# Distributed as part of the metahdep package, under the terms of the GNU General Public License (see DESCRIPTION file)

`metahdep.FEMA` <-
function(theta, V, X, meta.name="meta-analysis", center.X=FALSE)
{
  if (center.X)
    X <- center.columns(X)
  X <- metahdep.check.X(X)
  if (is.null(X))
    return(NULL)

  if( length(unique(c(length(theta),nrow(V),ncol(V),nrow(X)))) > 1 )
  {
    cat("Error:  Passed arguments are not conformable\n")
    return(NULL)
  }

  # Now run FEMA
  V.inv <- try(chol2inv(chol(V)))
  Sigma.beta <- try(chol2inv(chol(t(X) %*% V.inv %*% X )))

  if (!is.matrix(V.inv) | !is.matrix(Sigma.beta))
  {
    cat("Error:  The passed variance/covariance matrix (V) is not positive definite\n")
    return(NULL)
  }

  beta.hat <- Sigma.beta %*% t(X) %*% V.inv %*% theta
  Y.minus.X.beta <- theta - X %*% beta.hat

  ## Test of homogeneity / test of model mis-specification
  ## pp. 311 & 345 of Cooper & Hedges
  ## p. 172 of Hedges & Olkin
  Q <- t(Y.minus.X.beta) %*% V.inv %*% Y.minus.X.beta
  p.test <- 1 - pchisq(Q, (nrow(X)-ncol(X)))

  Z.beta <- abs(beta.hat/sqrt(diag(Sigma.beta)))
  P.beta <- 2*(1-pnorm(Z.beta))

  beta.hat <- t(beta.hat)
  dimnames(beta.hat) <- dimnames(X)
  P.beta <- t(P.beta)
  dimnames(P.beta) <- dimnames(X)

  return.list <- list(beta.hat, Sigma.beta, P.beta, Q, p.test, meta.name)
  names(return.list) <- c("beta.hats", "cov.matrix", "beta.hat.p.values", "Q", "Q.p.value", "name")

  return(return.list)
}


# Copyright 2008 John R. Stevens
# Distributed as part of the metahdep package, under the terms of the GNU General Public License (see DESCRIPTION file)

metahdep.list2dataframe <- function(output.list, max.rank, method, cov.names=NULL)
{
  num.genes <- length(output.list)

  ##  if length(cov.names) != max.rank then something is wrong.  use generic names instead
  if (length(cov.names) != max.rank)
    cov.names <- NULL

  ##  these will be the same no matter which method is used
  temp.bhatnames <- rep("", max.rank)
  temp.pvalnames <- rep("", max.rank)
  temp.covnames <- rep("", max.rank*max.rank)
  for (i in 1:max.rank)
  {
    if (!is.null(cov.names))
      temp.bhatnames[i] <- cov.names[i]
    else
      temp.bhatnames[i] <- paste("beta", (i-1), "hat", sep="")
    temp.pvalnames[i] <- paste(temp.bhatnames[i], ".pval", sep="")
    for (j in 1:max.rank)
      temp.covnames[j + (i-1)*max.rank] <- paste("cov", i, ".", j, sep="")
  }

  if (method == "FEMA")
  {
    ##  FEMA returns:  beta.hats, cov.matrix, beta.hat.p.values, Q, Q.p.value, gene
    name.mat <- rep("NA", num.genes)
    betahat.mat <- matrix(NA, nrow=num.genes, ncol=max.rank)
    covmat.mat <- matrix(NA, nrow=num.genes, ncol=max.rank*max.rank)
    betahat.pvals.mat <- matrix(NA, nrow=num.genes, ncol=max.rank)  
    q.mat <- matrix(NA, nrow=num.genes, ncol=1)
    qpval.mat <- matrix(NA, nrow=num.genes, ncol=1)
    for (i in 1:num.genes)
    {
      current.gene <- output.list[[i]]
      num.p <- length(current.gene[[1]])
      if (!is.matrix(current.gene[[2]]))
        current.gene[[2]] <- as.matrix(current.gene[[2]])
      if (num.p < max.rank)
      {
        temp.betahats <- rep(NA, max.rank)
        temp.betahat.pvals <- rep(NA, max.rank)
        temp.betahats[1:num.p] <- current.gene[[1]]
        temp.betahat.pvals[1:num.p] <- current.gene[[3]]

        temp.covmatrix <- matrix(NA, max.rank, max.rank)
        for (j in 1:num.p)
          temp.covmatrix[j,1:num.p] <- current.gene[[2]][j,]

        current.gene[[1]] <- temp.betahats
        current.gene[[2]] <- temp.covmatrix
        current.gene[[3]] <- temp.betahat.pvals
      }

      name.mat[i] <- current.gene[[6]]
      betahat.mat[i, 1:max.rank] <- current.gene[[1]]
      covmat.mat[i, 1:(max.rank*max.rank)] <- current.gene[[2]]
      betahat.pvals.mat[i, 1:max.rank] <- current.gene[[3]]
      q.mat[i] <- current.gene[[4]]
      if (is.na(q.mat[i]))
        qpval.mat[i] <- NA
      else
        qpval.mat[i] <- current.gene[[5]]
    }

    FEMA.frame <- cbind(betahat.mat, covmat.mat, betahat.pvals.mat, q.mat, qpval.mat)
    row.names(FEMA.frame) <- name.mat

    FEMA.colnames <- c(temp.bhatnames, temp.covnames, temp.pvalnames, "Q", "Q.p.value")
    FEMA.frame <- data.frame(FEMA.frame)
    names(FEMA.frame) <- FEMA.colnames
    return(FEMA.frame)
  }

  if (method == "REMA")
  {
    ##  REMA returns:  beta.hat, cov.matrix, P.values, sigma2, varsigma, Q, Q.p.value, gene
    name.mat <- rep("NA", num.genes)
    betahat.mat <- matrix(NA, nrow=num.genes, ncol=max.rank)
    covmat.mat <- matrix(NA, nrow=num.genes, ncol=max.rank*max.rank)
    betahat.pvals.mat <- matrix(NA, nrow=num.genes, ncol=max.rank)
    sigma2hat.mat <- matrix(NA, nrow=num.genes, ncol=1)
    varsigmahat.mat <- matrix(NA, nrow=num.genes, ncol=1)
    q.mat <- matrix(NA, nrow=num.genes, ncol=1)
    qpval.mat <- matrix(NA, nrow=num.genes, ncol=1)

    for (i in 1:num.genes)
    {
      current.gene <- output.list[[i]]
      num.p <- length(current.gene[[1]])
      if (!is.matrix(current.gene[[2]]))
        current.gene[[2]] <- as.matrix(current.gene[[2]])
      if (num.p < max.rank)
      {
        temp.betahats <- rep(NA, max.rank)
        temp.betahat.pvals <- rep(NA, max.rank)
        temp.betahats[1:num.p] <- current.gene[[1]]
        temp.betahat.pvals[1:num.p] <- current.gene[[3]]

        temp.covmatrix <- matrix(NA, max.rank, max.rank)
        for (j in 1:num.p)
          temp.covmatrix[j,1:num.p] <- current.gene[[2]][j,]

        current.gene[[1]] <- temp.betahats
        current.gene[[2]] <- temp.covmatrix
        current.gene[[3]] <- temp.betahat.pvals
      }

      name.mat[i] <- current.gene[[8]]
      betahat.mat[i, 1:max.rank] <- current.gene[[1]]
      covmat.mat[i, 1:(max.rank*max.rank)] <- current.gene[[2]]
      betahat.pvals.mat[i, 1:max.rank] <- current.gene[[3]]
      sigma2hat.mat[i,] <- current.gene[[4]]
      varsigmahat.mat[i,] <- current.gene[[5]]
      q.mat[i,] <- current.gene[[6]]
      if (is.na(q.mat[i,]))
        qpval.mat[i,] <- NA
      else
        qpval.mat[i,] <- current.gene[[7]]
    }

    REMA.frame <- cbind(betahat.mat, covmat.mat, betahat.pvals.mat, sigma2hat.mat, varsigmahat.mat, q.mat, qpval.mat)
    row.names(REMA.frame) <- name.mat
    REMA.colnames <- c(temp.bhatnames, temp.covnames, temp.pvalnames, "tau2.hat", "varsigma.hat", "Q", "Q.p.value")
    REMA.frame <- data.frame(REMA.frame)
    names(REMA.frame) <- REMA.colnames
    return(REMA.frame)
  }

  ##  else if (method == "HBLM")
  ##  HBLM returns:  betahats, var/covmat, betahat.pvals, sigma, sigma.var, varsigma, varsigma.var, sigma.varsigma.cov, gene
  name.mat <- rep("NA", num.genes)
  betahat.mat <- matrix(NA, nrow=num.genes, ncol=max.rank)
  covmat.mat <- matrix(NA, nrow=num.genes, ncol=max.rank*max.rank)
  betahat.pvals.mat <- matrix(NA, nrow=num.genes, ncol=max.rank)
  sigmahat.mat <- matrix(NA, nrow=num.genes, ncol=1)
  sigma.varhat.mat <- matrix(NA, nrow=num.genes, ncol=1)
  varsigmahat.mat <- matrix(NA, nrow=num.genes, ncol=1)
  varsigma.varhat.mat <- matrix(NA, nrow=num.genes, ncol=1)
  sigma.varsigma.cov.mat <- matrix(NA, nrow=num.genes, ncol=1)

  for (i in 1:num.genes)
  {
    current.gene <- output.list[[i]]
    num.p <- length(current.gene[[1]])
    if (!is.matrix(current.gene[[2]]))
      current.gene[[2]] <- as.matrix(current.gene[[2]])
    if (num.p < max.rank)
    {
      temp.betahats <- rep(NA, max.rank)
      temp.betahat.pvals <- rep(NA, max.rank)
      temp.betahats[1:num.p] <- current.gene[[1]]
      temp.betahat.pvals[1:num.p] <- current.gene[[3]]

      temp.covmatrix <- matrix(NA, max.rank, max.rank)
      for (j in 1:num.p)
        temp.covmatrix[j,1:num.p] <- current.gene[[2]][j,]

      current.gene[[1]] <- temp.betahats
      current.gene[[2]] <- temp.covmatrix
      current.gene[[3]] <- temp.betahat.pvals
    }

    name.mat[i] <- current.gene[[9]]
    betahat.mat[i, 1:max.rank] <- current.gene[[1]]
    covmat.mat[i, 1:(max.rank*max.rank)] <- current.gene[[2]]
    betahat.pvals.mat[i, 1:max.rank] <- current.gene[[3]]
    sigmahat.mat[i] <- current.gene[[4]]
    sigma.varhat.mat[i] <- current.gene[[5]]
    varsigmahat.mat[i] <- current.gene[[6]]
    varsigma.varhat.mat[i] <- current.gene[[7]]
    sigma.varsigma.cov.mat[i] <- current.gene[[8]]
  }

  HBLM.frame <- cbind(betahat.mat, covmat.mat, betahat.pvals.mat, sigmahat.mat, sigma.varhat.mat, varsigmahat.mat, varsigma.varhat.mat, sigma.varsigma.cov.mat)
  row.names(HBLM.frame) <- name.mat
  HBLM.colnames <- c(temp.bhatnames, temp.covnames, temp.pvalnames, "tau.hat", "tau.var", "varsigma.hat", "varsigma.var", "tau.varsigma.cov")
  HBLM.frame <- data.frame(HBLM.frame)
  names(HBLM.frame) <- HBLM.colnames
  return(HBLM.frame)

}

#####################################################################################
##  REMA Meta-analysis
#####################################################################################
LinMod.MetAn.dep.REMA <- function(gene.info)
{
  theta <- gene.info@theta
  V <- gene.info@V
  X <- gene.info@X
  M <- gene.info@M

  # Now run REMA
  temp.list <- LinMod.REMA.dep(gene.info)
  beta.hat <- temp.list[[1]]
  Sigma.beta <- temp.list[[2]]
  sig2.hat <- temp.list[[3]]
  varsig.hat <- temp.list[[4]]

  names(beta.hat) <- names(X)

  # REMA.testspec
  temp.W <- V + sig2.hat * id(nrow(V))
  resid <- theta - X %*% beta.hat
  Q <- as.numeric(t(resid) %*% chol2inv(chol(temp.W)) %*% resid)
  p.test <- 1 - pchisq(Q, (nrow(X)-ncol(X)))

  # REMA.testcovariates
  if(ncol(X)==1)
    t.beta <- abs(beta.hat/sqrt(Sigma.beta))
  else
    t.beta <- abs(beta.hat/sqrt(diag(Sigma.beta)))
  P.beta <- 2*(1-pt(abs(t.beta),df=nrow(X)-ncol(X)))

  beta.hat <- t(beta.hat)
  dimnames(beta.hat) <- dimnames(X)
  P.beta <- t(P.beta)
  dimnames(P.beta) <- dimnames(X)

  ##  make the return list
  return.list <- list(beta.hat, Sigma.beta, P.beta, sig2.hat, varsig.hat, Q, p.test, gene.info@gene)
  names(return.list) <- c("beta.hats", "cov.matrix", "beta.hat.p.values", "tau2.hat", "varsigma.hat", "Q", "Q.p.value", "name")

  return(return.list)
}     


###  used by LinMod.MetAn.dep.REMA to compute some of the parameter estimates
LinMod.REMA.dep <- function(gene.info)
{
  theta <- gene.info@theta
  V <- gene.info@V
  X <- gene.info@X
  M <- gene.info@M
  max.k <- gene.info@max.k

  if( length(unique(c(length(theta),nrow(V),ncol(V),nrow(X),nrow(M),nrow(M)))) > 1 )
  {
    cat("Error -- dimensions not compatible")
    return(NULL)
  }

  sig2.hat <- NA
  varsig.hat <- NA

  # The method of moments approach
  Vinv <- chol2inv(chol(V))
  inv.XVinvX <- chol2inv(chol(t(X) %*% Vinv %*% X))
  beta.hat <- inv.XVinvX %*% t(X) %*% Vinv %*% theta
  Y.minus.Xb <- theta - X %*% beta.hat
  RSS <- t(Y.minus.Xb) %*% Vinv %*% Y.minus.Xb

  sig2.hat <- as.numeric( (RSS-nrow(X)+ncol(X)) / sum(diag( Vinv - Vinv %*% X %*% inv.XVinvX %*% t(X) %*% Vinv )) )
  sig2.hat <- max(c(0,sig2.hat))
 
  psi <- V+sig2.hat*diag(c(rep(1,nrow(V))))

  inv.psi <- chol2inv(chol(psi))
  Sigma.b.w <- chol2inv(chol(t(X) %*% inv.psi %*% X))
  b.hat.w <- Sigma.b.w %*% t(X) %*% inv.psi %*% theta

  ##beta.Sigma.w <- cbind(b.hat.w,Sigma.b.w)
  beta.hat <- b.hat.w
  Sigma.beta <- Sigma.b.w

  return(list(as.vector(beta.hat), Sigma.beta, sig2.hat, varsig.hat))
}


##################################################################################
### new version of REMA with delta splitting
### Function to perform Random Effects meta-analysis delta-splitting
### Arguments: 
### gene.info = prep.list element containing the collected information for the current gene
###             this contains:  theta, V, X, M, max.k, and the gene name
### epsilon = convergence criterion for tausq(sigma2) and varsigma
### maxiter = maximum number of iterations for REMA estimation
##################################################################################
LinMod.REMA.delta.split <- function(gene.info, epsilon=1e-05, maxiter=100)
{
  theta <- gene.info@theta
  V <- gene.info@V
  X <- gene.info@X
  M <- gene.info@M
  K <- gene.info@max.k

  num <- maxiter # number of iterations
  sigma2.vec <- varsigma.vec <- RSS.sigma2.vec <- RSS.varsigma.vec <- rep(NA,num)
  # 1.
  sigma2 <- varsigma <- 0

  keep.iterating <- TRUE
  i <- 0
  Y <- theta
  while(keep.iterating)
  {
    i <- i+1
    # 2. Get A and RSS
    psi.inv <- chol2inv(chol(V + sigma2*id(nrow(V)) + varsigma*M))
    A <- psi.inv - psi.inv %*% X %*% chol2inv(chol(t(X) %*% psi.inv %*% X)) %*% t(X) %*% psi.inv
    RSS <- as.numeric(t(Y) %*% A %*% Y)
    RSS.sigma2.vec[i] <- RSS
    # 3. Update tau (sigma2)
    sigma2 <- as.numeric((RSS - tr(A %*% V) - varsigma*tr(A %*% M)) / tr(A))
    if(sigma2<0)
      {  sigma2 <- 0  }
    sigma2.vec[i] <- sigma2
    if (varsigma < -sigma2/(K - 1)) { varsigma <- (-sigma2/(K - 1)) }   # added 11.25.08
    if (varsigma > sigma2)  { varsigma <- sigma2 }                      # added 11.25.08
    # 4. Get Psi
    psi.inv <- chol2inv(chol((V + sigma2*id(nrow(V)) + varsigma*M)))
    # 5. Get A and RSS
    A <- psi.inv - psi.inv %*% X %*% chol2inv(chol(t(X) %*% psi.inv %*% X)) %*% t(X) %*% psi.inv

    RSS <- as.numeric(t(Y) %*% A %*% Y)
    RSS.varsigma.vec[i] <- RSS
    # 6. Update varsigma
    varsigma <- (RSS - tr(A %*% V) - sigma2*tr(A)) / tr(A %*% M)

    ##  check for valid results; if there's a problem then stop iterating
    if (is.na(sigma2) | is.nan(sigma2) | is.na(varsigma) | is.nan(varsigma))
    {
      if (is.na(sigma2) | is.nan(sigma2))
        sigma2 <- NA
      varsigma <- NA
      keep.iterating <- FALSE
    }    else
    {
      if(varsigma < -sigma2/(K-1))
        {  varsigma <-  (-sigma2/(K-1))  }
      if(varsigma > sigma2)
        {  varsigma <- sigma2  }
      varsigma.vec[i] <- varsigma
      # Check convergence beginning at fifth iteration
      if(i > 4)
        if( abs(sigma2.vec[i]-sigma2.vec[i-1]) < epsilon & abs(varsigma.vec[i]-varsigma.vec[i-1]) < epsilon)
          keep.iterating <- FALSE
      if(i >= num)
        keep.iterating <- FALSE
    }
  }

  beta.hat <- chol2inv(chol(t(X) %*% psi.inv %*% X)) %*% t(X) %*% psi.inv %*% Y
  Sigma.beta <- chol2inv(chol(t(X) %*% psi.inv %*% X))
  se.beta.hat <- sqrt(diag(Sigma.beta))
  t.beta <- abs(beta.hat/se.beta.hat)
  P.beta <- 2*(1-pt(abs(t.beta),df=nrow(X)-ncol(X)))

  ##  compute the Q statistic and it's p-value
  resid <- theta - X %*% beta.hat
  Q <- as.numeric(t(resid) %*% psi.inv %*% resid)
  p.test <- 1 - pchisq(Q, (nrow(X)-ncol(X)))

  beta.hat <- t(beta.hat)
  dimnames(beta.hat) <- dimnames(X)
  P.beta <- t(P.beta)
  dimnames(P.beta) <- dimnames(X)

  return.list <- list(beta.hat, Sigma.beta, P.beta, sigma2, varsigma, Q, p.test, gene.info@gene)
  names(return.list) <- c("beta.hats", "cov.matrix", "beta.hat.p.values", "tau2.hat", "varsigma.hat", "Q", "Q.p.value", "name")
  return(return.list)
}



#####################################################################################
##  HBLM Meta-analysis
#####################################################################################
### Fit an HBLM (without delta-splitting)
### This returns a five-part list:
###    [[1]] is the vector of covariate estimates
###    [[2]] is the corresponding covariance matrix
###    [[3]] is the corresponding vector of significance P-values,
###    [[4]] is the posterior mean of the variance sigma
###    [[5]] is the posteror variance of the variance sigma.
LinMod.HBLM.fast.dep <- function(gene.info, n=10, two.sided=FALSE)
{
  theta <- gene.info@theta
  V <- gene.info@V
  X <- gene.info@X

  # Make theta the right kind of object
  theta <- matrix(theta,ncol=1)

  # Get harmonic mean of sampling variances
  c.0 <- sqrt(nrow(V)/sum(1/diag(V)))

  # Define steps to take in integrations; make steps and step.size global variables.
  first.steps <- c(0:n)*(c.0/3)/n
  second.steps <- (c.0/3) + c(0:n)*(2*c.0/3)/n
  third.steps <- (c.0)+c(0:n)*(2*c.0/n)
  fourth.steps <- 3*c.0 + c(0:n)*c.0
  steps <- matrix(c(first.steps,second.steps,third.steps,fourth.steps),ncol=1)
  base.step <- (2^(c(0,rep(1,(n-1)),0)+mod(c(0:n),2)))/3   # This generates c(1,4,2,4,2,4,...,2,4,1)/3
  step.size <- c.0*c(base.step/3/n,2*base.step/3/n,2*base.step/n,base.step)
  sigma.v <- steps
  delta.v <- step.size

  # Get array of matrices Psi = (V+sigma^2*I)^{-1}
  Psi.v <- array(dim=c(nrow(V),ncol(V),length(sigma.v)))
  I.v <- id(nrow(V))
  for(i in 1:length(sigma.v))
    Psi.v[,,i] <- chol2inv(chol(V+sigma.v[i]^2*I.v))

  # Get array of matrices beta.star.star = (X^T Psi X)^{-1}
  beta.star.star.v <- array(dim=c(ncol(X),ncol(X),dim(Psi.v)[3]))
  for(i in 1:(dim(Psi.v)[3]))
    beta.star.star.v[,,i] <- chol2inv(chol(t(X) %*% Psi.v[,,i] %*% X))

  # Get array of matrices S = Psi - Psi X beta.star.star t(X) Psi
  S.v <- array(dim=dim(Psi.v))
  for(i in 1:(dim(Psi.v)[3]))
    S.v[,,i] <- Psi.v[,,i] - Psi.v[,,i] %*% X %*% beta.star.star.v[,,i] %*% t(X) %*% Psi.v[,,i]

  # Get matrix with columns beta.star = beta.star.star t(X) Psi theta
  beta.star.v <- matrix(nrow=ncol(X),ncol=dim(Psi.v)[3])
  for(i in 1:(dim(Psi.v)[3]))
    beta.star.v[,i] <- beta.star.star.v[,,i] %*% t(X) %*% Psi.v[,,i] %*% theta

  # Get vector of prior values pi = c/(c+sigma)^2
  pi.prior.v <- c.0/(c.0+sigma.v)^2

  # Get vector of values proportional to posterior : p  ~[proportional to]~ pi.posterior(sigma|theta)
  p.v <- rep(NA,length(pi.prior.v))
  for(i in 1:length(pi.prior.v))
    p.v[i] <- pi.prior.v[i] * prod(diag(chol(Psi.v[,,i]))) * prod(diag(chol(as.matrix(beta.star.star.v[,,i])))) * exp(-.5*as.numeric(t(theta) %*% S.v[,,i] %*% theta))

  # Get normalizing scalar value for posterior
  norm.value <- as.numeric(t(p.v) %*% delta.v)

  # Get posterior mean of beta given theta
  beta.star <- as.vector(1/norm.value * beta.star.v %*% (p.v*delta.v))

  # Get posterior covariance of beta given theta
  B.v <- array(dim=dim(beta.star.star.v))
  for(i in 1:(dim(B.v)[3]))
    B.v[,,i] <- ( beta.star.star.v[,,i] + (beta.star.v[,i]-beta.star) %*% t(beta.star.v[,i]-beta.star) ) * p.v[i] * delta.v[i]
  beta.star.star <- 1/norm.value * apply(B.v,c(1,2),sum)

  # Get phi values for use in computing posterior probabilities
  phi.v <- matrix(nrow=nrow(beta.star.v),ncol=ncol(beta.star.v))
  if(length(beta.star.star.v[,,1])==1)
    phi.v <- matrix( pnorm(as.vector(beta.star.v)/as.vector(beta.star.star.v)) ,nrow=1)
  else
  {
    for(i in 1:ncol(phi.v))
      phi.v[,i] <- pnorm(beta.star.v[,i] / sqrt(diag(beta.star.star.v[,,i]))) 
  }

  # Get posterior probabilities P(beta[j] > 0 | data)
  beta.pval <- as.vector(1/norm.value * phi.v %*% (p.v*delta.v))

  # if two.sided==TRUE then get P(beta[j] != 0 | data) instead
  if (two.sided)
    beta.pval <- 1-2*abs(.5-beta.pval)

  # Get estimated posterior mean of variance : E[sigma | data]
  E.post.sigma <- as.numeric(1/norm.value * t(sigma.v) %*% (p.v*delta.v))

  # Get estimated posterior variance of variance : Var[sigma | data]
  E.post.sigma.squared <- as.numeric(1/norm.value * t(sigma.v^2) %*% (p.v*delta.v))
  V.post.sigma <- E.post.sigma.squared - E.post.sigma^2

  # Return results as a list
  return.list <- list(beta.star, beta.star.star, beta.pval, E.post.sigma, V.post.sigma, NA, NA, NA, gene.info@gene)
  names(return.list[[1]]) <- colnames(X)
  names(return.list) <- c("beta.hats", "cov.matrix", "beta.hat.p.values", "tau.hat", "tau.var", "varsigma.hat",
                              "varsigma.var", "tau.varsigma.cov", "name")
  return(return.list)
}



##  fit an HBLM with delta-splitting by using numerical integration over sigma and varsigma with Simpsons rule
##  the integration over sigma is actually four applications of Simpsons rule -- done for each quartile
##
##  the new function has also been updated to properly account for dependence groups (via the matrix gene.m)
##  and also to use a new lower bound when integrating varsigma (-sigma^2/(k-1) instead of zero)
new.LinMod.HBLM.fast.dep.delta.split <- function(gene.info, n=10, m=10, two.sided=FALSE)
{
    theta <- gene.info@theta
    V <- gene.info@V
    X <- gene.info@X
    M <- gene.info@M
    max.k <- gene.info@max.k

    # Make theta the right kind of object
    theta <- matrix(theta,ncol=1)
    if (!is.matrix(X))
      X <- as.matrix(X)
    tX <- t(X)
    I.v <- diag(rep(1,nrow(V)))

    # check max.k
    if (as.integer(max.k) < 2)
    {
      cat("Error:  max.k must be an integer greater than 1 -- current value is", max.k, "\n")
      return(NULL)
    }
    max.k <- as.integer(max.k)

    # Get harmonic mean of sampling variances
    c.0 <- sqrt(nrow(V)/sum(1/diag(V)))

    first.steps <- c(0:n)*(c.0/3)/n
    second.steps <- (c.0/3) + c(0:n)*(2*c.0/3)/n
    third.steps <- (c.0)+c(0:n)*(2*c.0/n)
    fourth.steps <- 3*c.0 + c(0:n)*c.0
    steps <- matrix(c(first.steps,second.steps,third.steps,fourth.steps),ncol=1)
    q.n <- (2^(c(0,rep(1,(n-1)),0)+mod(c(0:n),2)))  # typo fixed 10.20.08 (last 2 was n)


    # Get necessary components (as in 12.16.04 notes)
    sigma.v <- steps
    delta.sigma <- c.0*c(rep(1/3/n,n+1),rep(2/3/n,n+1),rep(2/n,n+1),rep(1,n+1))
    varsigma.v <- t(apply(sigma.v,1,get.varsigma.v,m=m,max.k=max.k))
    delta.varsigma <- 2*sigma.v/m
    Q.Delta <- (rep(q.n,4)*delta.varsigma/3*delta.sigma/3)%*%t((2^(c(0,rep(1,(m-1)),0)+mod(c(0:m),2))))

    # Get vector of prior values pi = 1/sigma^2 * c/(c+sigma)^2
    pi.prior.v <- 1/(.0000001+sigma.v^2) * c.0/(c.0+sigma.v)^2


    num.sigma.steps <- nrow(varsigma.v)      ##  the number of abscissas for evaluating sigma
    num.varsigma.steps <- ncol(varsigma.v)      ##  the number of abscissas for evaluating varsigma

    Qp <- array(0., dim=c(num.sigma.steps, num.varsigma.steps))
    beta.star.star.v <- array(dim=c(num.sigma.steps, num.varsigma.steps, ncol(X), ncol(X)))
    beta.star.v <- array(dim=c(num.sigma.steps, num.varsigma.steps, ncol(X)))
    J.twiddle.star <- array(dim=c(num.sigma.steps, num.varsigma.steps, ncol(X)))


    ##  step through the grid of points, evaluating the function value at each point
    ##  save the appropriate data
    for(i in 1:num.sigma.steps)
    {
      for(j in 1:num.varsigma.steps)
      {
        ##  Psi = (V+sigma^2*I+varsigma*M)^{-1}
        Psi <- chol2inv(chol(V + sigma.v[i]^2*I.v + varsigma.v[i,j]*M))

        ##  beta.star.star = (X^T Psi X)^{-1}
        ##  save this matrix, as it is used later
        beta.star.star <- chol2inv(chol((tX %*% Psi %*% X)))
        beta.star.star.v[i,j,,] <- beta.star.star

        ##  S = Psi - Psi X beta.star.star t(X) Psi
        S <- Psi - Psi %*% X %*% beta.star.star %*% tX %*% Psi

        ##  Get matrix with columns beta.star = beta.star.star t(X) Psi theta
        ##  save this matrix, as it is used later
        beta.star.v[i,j,] <- beta.star.star %*% tX %*% Psi %*% theta

        ##  Get value proportional to posterior : p  ~[proportional to]~ pi.posterior(sigma,varsigma|theta)
        ##p <- pi.prior.v[i] * sqrt(det(Psi)) * sqrt(det(as.matrix(beta.star.star))) * exp(-.5*as.numeric(t(theta) %*% S %*% theta))
        p <- pi.prior.v[i] * prod(diag(chol(Psi))) * prod(diag(chol(beta.star.star))) * exp(-.5*as.numeric(t(theta) %*% S %*% theta))

        ##  update a vector that is used to compute the normalizing value
        Qp[i,j] <- Q.Delta[i,j] * p
      }

      ##  this will also be used to compute the normalizing value
      J.twiddle.star[i,,] <- Qp[i,]
    }

    # Get normalizing scalar value for posterior
    norm.value <- sum(Qp)


    ##  get the vector of coefficient estimates (beta-hats)
    J.b <- J.twiddle.star*beta.star.v
    beta.star <- (1/norm.value) * apply(J.b,3,sum)


    # Get posterior covariance of beta given theta
    J.twiddle.star.star <- array(dim=c(num.sigma.steps, num.varsigma.steps, ncol(X), ncol(X)))
    B.diff <- array(dim=c(num.sigma.steps, num.varsigma.steps, ncol(X), ncol(X)))
    for(i in 1:length(sigma.v))
      for(j in 1:ncol(varsigma.v))
      {
        J.twiddle.star.star[i,j,,] <- Qp[i,j]
        B.diff[i,j,,] <- (beta.star.v[i,j,] - beta.star) %*% t(beta.star.v[i,j,] - beta.star)
      }
    J.B <- J.twiddle.star.star*(beta.star.star.v+B.diff)
    beta.star.star <- (1/norm.value) * apply(J.B,c(3,4),sum)


    # Get phi values for use in computing posterior probabilities
    phi.v <- array(dim=dim(beta.star.v))
    if(ncol(X)==1)
      for(i in 1:num.sigma.steps)
        for(j in 1:num.varsigma.steps)
          phi.v[i,j,] <- pnorm(as.numeric(beta.star.v[i,j,]) / sqrt( as.numeric((beta.star.star.v[i,j,,]))))
    else
      for(i in 1:num.sigma.steps)
        for(j in 1:num.varsigma.steps)
          phi.v[i,j,] <- pnorm(beta.star.v[i,j,] / sqrt(diag(beta.star.star.v[i,j,,])))


    # Get posterior probabilities P(beta[j] > 0 | data)
    J.PHI <- J.twiddle.star*phi.v
    beta.pval <- (1/norm.value) * apply(J.PHI, 3, sum)

    # if two.sided==TRUE then get P(beta[j] != 0 | data) instead
    if (two.sided)
      beta.pval <- 1-2*abs(.5-beta.pval)


    # Get estimated posterior mean of variance : E[sigma | data]
    sigma.matrix <- (matrix(rep(sigma.v, (m+1)), ncol=(m+1)))
    E.post.sigma <- (1/norm.value) * sum(Qp*sigma.matrix)


    # Get estimated posterior variance of variance : Var[sigma | data]
    E.post.sigma.squared <- (1/norm.value) * sum(Qp*sigma.matrix*sigma.matrix)
    V.post.sigma <- E.post.sigma.squared - E.post.sigma^2


    # Get estimated posterior mean and variance of covariance
    E.post.varsigma <- (1/norm.value) * sum(Qp*varsigma.v)
    E.post.varsigma.squared <- (1/norm.value) * sum(Qp*varsigma.v*varsigma.v)
    V.post.varsigma <- E.post.varsigma.squared - E.post.varsigma^2


    # Get estimated posterior covariance of the variance and covariance
    E.post.sigma.varsigma <- (1/norm.value) * sum(Qp*sigma.matrix*varsigma.v)
    Cov.post.sigma.varsigma <- E.post.sigma.varsigma - E.post.sigma * E.post.varsigma


    # Return results as a list
    return.list <- list(beta.star, beta.star.star, beta.pval, E.post.sigma, V.post.sigma, E.post.varsigma,
                             V.post.varsigma, Cov.post.sigma.varsigma, gene.info@gene)
    names(return.list[[1]]) <- colnames(X)
    names(return.list) <- c("beta.hats", "cov.matrix", "beta.hat.p.values", "tau.hat", "tau.var", "varsigma.hat",
                              "varsigma.var", "tau.varsigma.cov", "name")
    return(return.list)
}


#####################################################################################
##  FEMA Meta-analysis
#####################################################################################

LinMod.MetAn.dep.FEMA <- function(gene.info)
{
  theta <- gene.info@theta
  V <- gene.info@V
  X <- gene.info@X
  max.k <- gene.info@max.k

  if( length(unique(c(length(theta),nrow(V),ncol(V),nrow(X)))) > 1 )
    return(NULL)

  # Now run FEMA
  V.inv <- try(chol2inv(chol(V)))
  Sigma.beta <- try(chol2inv(chol(t(X) %*% V.inv %*% X )))

  if (!is.matrix(V.inv) | !is.matrix(Sigma.beta))
    return(NULL)

  beta.hat <- Sigma.beta %*% t(X) %*% V.inv %*% theta

  Y.minus.X.beta <- theta - X %*% beta.hat

  ## Test of homogeneity / test of model mis-specification
  ## pp. 311 & 345 of Cooper & Hedges
  ## p. 172 of Hedges & Olkin
  Q <- t(Y.minus.X.beta) %*% V.inv %*% Y.minus.X.beta
  p.test <- 1 - pchisq(Q, (nrow(X)-ncol(X)))

  t.beta <- abs(beta.hat/sqrt(diag(Sigma.beta)))
  P.beta <- 2*(1-pt(abs(t.beta),df=nrow(X)-ncol(X)))

  beta.hat <- t(beta.hat)
  dimnames(beta.hat) <- dimnames(X)
  P.beta <- t(P.beta)
  dimnames(P.beta) <- dimnames(X)

  return.list <- list(beta.hat, Sigma.beta, P.beta, Q, p.test, gene.info@gene)
  names(return.list) <- c("beta.hats", "cov.matrix", "beta.hat.p.values", "Q", "Q.p.value", "gene")

  return(return.list)
}     



#####################################################################################
##  Miscellaneous functions
#####################################################################################

#######  check to see if a matrix X has too many columns
#######  i.e. if there are N observations in theta, and P covariates in X, then X will be NxP
#######  if P >= N then the there are too many things being estimated and not enough observations
#######  so, this function cuts columns from X until P < N
#######  also, it check that the columns are linearly independent
metahdep.check.X <- function(X)
{
  if (!is.matrix(X))
    X <- as.matrix(X)
  N <- nrow(X)
  P <- ncol(X)
  if (N < 2)
  {
    cat("Error:  Cannot perform a meta-analysis on a single observation (nrow(X) < 2)\n")
    return(NULL)
  }
  X.rank <- qr(X)$rank
  if (X.rank < P)
  {
    cat("Error:  The passed covariate matrix (X) is not full rank\n")
    return(NULL)
  }
  if (P == N)
  {
    temp.X <- X[,1:(P-1)]
    X <- temp.X
    cat("Warning:  The covariate matrix (X) had too many columns -- a covariate has been removed\n")
  }
  return(X)
}

####### Obtain M matrix, given information about the dependence group structure
##  dep.groups is a list of vectors of study ids, like:
##  dep.groups <- list(c(2,3,4,5),c(9,10,11))
get.M <- function(N, dep.groups)
{
  if (is.matrix(N))
  {
    M <- 0*N
    N <- nrow(M)
  }
  else
    M <- matrix(0, N, N)

  ##  check dep.groups
  N.temp <- 0
  for (i in 1:length(dep.groups))
    N.temp <- N.temp + length(dep.groups[[i]])
  if (N.temp > N)
  {
    cat("Error:  the passed dep.groups object indicates a dependence matrix M that is too big\n")
    cat("Expected size: ", N, "x", N, "\n")
    cat("Size indicated by dep.groups: ", N.temp, "x", N.temp, "\n")
    return(NULL)
  }

  k.vec <- rep(NA,length(dep.groups))
  for(i in 1:length(dep.groups))
  {
    dep.studies <- dep.groups[[i]]
    k.vec[i] <- length(dep.studies)
    use.grid <- as.matrix(expand.grid(dep.studies,dep.studies))

    #throw out diagonal values
    use.grid <- use.grid[!(use.grid[,1]==use.grid[,2]),]
    use.grid <- use.grid[(!(use.grid[,1] > N) & !(use.grid[,2] > N)),]
    M[use.grid] <- 1
  }
  return(M)
}


# get the trace of a matrix
tr <- function(B)
{
  if (is.matrix(B))
    return(sum(diag(B)))
  else
    return(sum(diag(as.matrix(B))))
}

# Create an identity matrix of dimension p
id <- function(p)
{  return(diag(rep(1,p)))  }

# Given a matrix X, this will return a matrix whose jth column is X[,j]-mean(X[,j]) for j > 1
# This is used in the metahdep functions (FEMA, REMA, HBLM) -- to make the 
# interpretation of the intercept term be the population effect size.
center.columns <- function(X)
{
  if (!is.matrix(X))
    X <- as.matrix(X)
  X.new <- X
  ncol.X <- ncol(X)
  if(ncol.X > 1)
    for(j in 2:ncol.X)
      X.new[,j] <- X[,j] - mean(X[,j])
  return(X.new)
}


##  Define mod function
##  This is redundant.  It can be replaced with the standard R operator:  %%
##  i.e. this should be equivalent to:  mod <- function(a,b) return(a %% b)
mod <- function(a,b)
{  return(a-b*trunc(a/b))  }


##  This is used by the fast HBLM delta splitting function via apply()
##  return values from -x^2/(max.k-1) to x^2 in m steps
get.varsigma.v <- function(x,m,max.k=4)
{
  lower.val <- -x^2/(max.k-1)
  upper.val <- x^2
  ret.seq <- seq(from=lower.val,to=upper.val,length.out=m+1)
  return(ret.seq)
}

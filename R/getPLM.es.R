# Copyright 2008 John R. Stevens
# Distributed as part of the metahdep package, under the terms of the GNU General Public License (see DESCRIPTION file)

`getPLM.es` <- function(abatch,trt1,trt2,covariates=NULL,dep.grp=NULL,sub.gn=NULL,bg.norm=TRUE)
  {
     require(affyPLM)
     # check arguments
     if(class(abatch)!="AffyBatch"){cat('Error: abatch must be AffyBatch','\n'); return(NULL) }
     if(class(trt1)!=class(trt2) ){cat('Error: trt1 and trt2 must be both lists or both numeric vectors','\n'); return(NULL)}
     if(class(trt1)=='list' & length(trt1)!=length(trt2)){cat('Error: if lists, trt1 and trt2 must be of same length','\n'); return(NULL)}
     if(!is.list(trt1)){ trt1 <- list(trt1); trt2 <- list(trt2) }
     if(mean(is.element(c(unlist(trt1),unlist(trt2)),1:length(abatch)))!=1){cat('Error: invalid array indices','\n'); return(NULL) }
     if(!is.null(covariates))
       { if(!is.data.frame(covariates))
            {cat('Error: covariates argument must be a data.frame object','\n'); return(NULL) } 
         if(is.data.frame(covariates) & nrow(covariates)!=length(trt1) ){cat('Error: number of rows in covariates argument must equal number of comparisons represented in trt1 and trt2 arguments','\n'); return(NULL)}
         if(is.data.frame(covariates)){ if(!is.numeric(as.matrix(covariates))){cat('Error: covariates values must be numeric','\n'); return(NULL)} }
        }
     if(!is.null(dep.grp)){ if(!is.numeric(dep.grp) | length(dep.grp)>1 | dep.grp!=as.integer(dep.grp) ){cat('Error: dep.grp argument must be a single integer value','\n'); return(NULL) } }

     # run fitPLM
     set <- fitPLM(abatch,PM~-1+samples+probes,background=bg.norm,normalize=bg.norm,subset=sub.gn,output.param=list(varcov="chiplevel"))
     eset <- coefs(set)
     vc <- varcov(set)

     # Get effect size estimates and covariances for comparisons of interest
     num.comp <- length(trt1)
     num.gn <- nrow(eset)
     num.covar <- num.comp*(num.comp+1)/2
     es <- matrix(nrow=num.gn,ncol=num.comp)
     v <- NULL
     for(j in 1:num.comp)
        {
          n1j <- length(trt1[[j]])
          n2j <- length(trt2[[j]])
          cj <- matrix(0,nrow=1,ncol=ncol(eset))
          cj[1,trt1[[j]]] <- (-1/n1j)
          cj[1,trt2[[j]]] <-  (1/n2j)
          fcj <- eset%*%t(cj)
          t.plm <- rep(0,num.gn)
          for(i in 1:num.gn)
             {
               t.plm[i] <- fcj[i]/sqrt(cj%*%vc[[i]]%*%t(cj))
              }
          es[,j] <- sqrt(1/n1j+1/n2j)*t.plm
          for(h in j:num.comp)
             {
               n1h <- length(trt1[[h]])
               n2h <- length(trt2[[h]])
               ch <- matrix(0,nrow=1,ncol=ncol(eset))
               ch[1,trt1[[h]]] <- (-1/n1h)
               ch[1,trt2[[h]]] <-  (1/n2h)
               v.plm <- denom.plm <- rep(0,num.gn)
               for(i in 1:num.gn)
                  {
                    v.plm[i] <- cj%*%vc[[i]]%*%t(ch)
                    denom.plm[i] <- cj%*%vc[[i]]%*%t(cj) * ch%*%vc[[i]]%*%t(ch)
                   }
               v.k <- sqrt((1/n1j+1/n2j)*(1/n1h+1/n2h)/denom.plm) * v.plm
               v <- cbind(v,v.k)               
              }
        }

     ## Put it all together in one ES.obj object and return this object
     gn <- as.character(rownames(eset))
     #
     es.names <- rep(NA,num.comp)
     for(i in 1:num.comp){es.names[i] <- paste('ES',(i-1),sep='.')}
     colnames(es) <- es.names
     #
     v.names <- NULL
     for(j in 1:num.comp)
        {
           for(h in j:num.comp)
              {
                v.names <- c(v.names,paste('Cov',(j-1),(h-1),sep='.'))
               }
         }
     colnames(v) <- v.names
     #
     result <- new("ES.obj", gn=gn, ES.mat=es, Cov.mat=v, chip=annotation(abatch) )
     if(!is.null(covariates)){result@covariates <- covariates }
     if(!is.null(dep.grp)){result@dep.grp <- as.integer(dep.grp) }
     return(result)
   }

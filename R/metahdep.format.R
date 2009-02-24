# Copyright 2008 John R. Stevens
# Distributed as part of the metahdep package, under the terms of the GNU General Public License (see DESCRIPTION file)

`metahdep.format` <-
function(ES.obj.list, newnames, min.var = 0.0001, 
                include.row.indices = FALSE, show.warnings = FALSE, 
                pd.verify = FALSE)
{
  ###########################################################
  ##	INITIALIZATION
  ###########################################################

  ## Make initial objects and check for appropriate formatting
  if(!is.list(ES.obj.list)){ cat('Error: ES.obj.list argument must be a list of ES.obj objects','\n'); return(NULL) }
  if(mean(unlist(lapply(ES.obj.list,class))=="ES.obj")!=1){ cat('Error: elements of ES.obj.list argument must be of class ES.obj','\n'); return(NULL) }

  # make dependence.groups object
  num.studies <- length(ES.obj.list)
  dependence.groups <- c()
  full.chipset.list <- c()
  for(i in 1:num.studies)
     {
        dependence.groups <- c(dependence.groups, ES.obj.list[[i]]@dep.grp)
        full.chipset.list <- c(full.chipset.list, ES.obj.list[[i]]@chip)
      }
  if(length(dependence.groups)<num.studies){ dependence.groups <- 1:num.studies; cat('Caution: dependence.groups information missing for at least one study ...','\n','... proceeding with dependence.groups defined as ',dependence.groups,'\n') }

  # make covariates object
  num.studies <- length(ES.obj.list)
  covar.names <- NULL
  num.comp <- 0
  for(i in 1:num.studies)
    {
       covar <- ES.obj.list[[i]]@covariates
       covar.names <- union(covar.names, names(covar))
       num.comp <- num.comp + nrow(covar)
     }
  Intercept <- rep(1,num.comp)
  design.matrix <- matrix(Intercept,ncol=1)
  for(i in 1:length(covar.names))
    {
       c.name <- covar.names[i]
       c.vector <- c()
       for(j in 1:num.studies)
          {
             covar <- ES.obj.list[[j]]@covariates
             if(!is.element(c.name,names(covar))){cat('Error: covariate ',c.name,' not defined in all studies','\n'); return(NULL) }
             t <- names(covar)==c.name
             c.vector <- c(c.vector,covar[,t])
           }
       design.matrix <- cbind(design.matrix,c.vector)     
     }
  colnames(design.matrix) <- c('Intercept',covar.names)

  covariates <- as.matrix(design.matrix)
  covariate.names <- colnames(covariates)
  p <- length(covariate.names)
  num.dep.groups <- length(unique(dependence.groups))
  pd.check <- as.matrix(1)

  ##  find the number of differential expression measures in each study
  ##  this may change if the studies are reordered, and so this is also done later in the function
  num.DifExpMeasures <- c(rep(0,num.studies))
  for(i in 1:num.studies)
  {
    num.DifExpMeasures[i] <- ncol(ES.obj.list[[i]]@ES.mat)
  }

  if (length(covariates[,1]) != sum(num.DifExpMeasures))
  {
    cat("Error:  The covariate matrix has an incorrect number of rows\n")
    cat("Rows expected: ", sum(num.DifExpMeasures), "\n")
    cat("Rows found: ", length(covariates[,1]), "\n")
    return(NULL)
  }

  ##  check the dependence.groups vector to make sure it is valid
  dependence.groups <- as.integer(dependence.groups)
  if (sum(is.na(dependence.groups)) > 0)
  {
    cat("Error:  dependence.groups argument must be a vector of integers or coercible into integers\n")
    return(NULL)
  }



  ##  check whether gn, ES.mat, and Cov.mat objects have right dimensions
  for(i in 1:num.studies)
    {
       num.comp <- ncol(ES.obj.list[[i]]@ES.mat)
       length.upper.Triangle <- num.comp*(num.comp+1)/2
       if( ncol(ES.obj.list[[i]]@Cov.mat) != length.upper.Triangle ){ cat('Error: dimensions of Cov.mat is incorrect for number of comparisons defined in study ',i,' of ES.obj.list argument','\n'); return(NULL) }
       if( sd( c( nrow(ES.obj.list[[i]]@ES.mat), nrow(ES.obj.list[[i]]@Cov.mat), length(ES.obj.list[[i]]@gn) ) ) != 0 ){cat('Error: ES.mat, Cov.mat, and gn slots in ES.obj do not have same number of rows in study ',i,' of ES.obj.list argument','\n'); return(NULL) }
     }

  ## Make the DifExp.list object - a list of data.frame objects, where each data.frame represents an ES.obj object
  DifExp.list <- list()
  for(i in 1:num.studies)
     {
        gn <- ES.obj.list[[i]]@gn
        num.genes <- length(gn)
        num.comp <- ncol(ES.obj.list[[i]]@ES.mat)
        label <- c(1:num.comp,rep('.',num.genes-num.comp))
        d.frame <- data.frame(label=label, gn=gn)
        DifExp.list[[i]] <- as.data.frame(cbind(d.frame,ES.obj.list[[i]]@ES.mat,ES.obj.list[[i]]@Cov.mat))
      }

  ##  need to rearrange the data frames in the DifExp.list object to be in the same order as sorted dependence.groups, 
  ##  as the function relies on this being the case when it contructs the return objects (M in particular)
  ##  the chipset name list and the rows of the covariate matrix will also need to be reordered
  dep.group.mat <- cbind(dependence.groups, (1:num.studies))
  dep.group.sort.indices <- sort(dependence.groups, index.return=TRUE)
  sorted.dep.group.mat <- dep.group.mat[dep.group.sort.indices$ix,]
  dependence.groups <- sorted.dep.group.mat[,1]
  temp.difexp.list <- list()
  for(i in 1:num.studies)
    temp.difexp.list[[i]] <- DifExp.list[[sorted.dep.group.mat[i,2]]]
  DifExp.list <- temp.difexp.list
  full.chipset.list <- full.chipset.list[dep.group.sort.indices$ix]

  # reorder the rows of the covariate matrix so they line up with the sorted DifExp.list
  cov.mat.offsets <- rep(0, num.studies)
  current.study.list <- dep.group.sort.indices$ix
  cov.mat.row.indices <- rep(0, sum(num.DifExpMeasures))
  cov.mat.row.indices.offset <- 0
  for (i in 2:num.studies)
    cov.mat.offsets[i] <- sum(num.DifExpMeasures[1:(i-1)])
  for (i in 1:num.studies)
  {
    current.num.dif.measures <- num.DifExpMeasures[current.study.list[i]]
    for (j in 1:current.num.dif.measures)
      cov.mat.row.indices[cov.mat.row.indices.offset + j] <- cov.mat.offsets[current.study.list[i]] + j
    cov.mat.row.indices.offset <- cov.mat.row.indices.offset + num.DifExpMeasures[current.study.list[i]]
  }
  temp.covariates <- covariates[cov.mat.row.indices,]
  covariates <- temp.covariates

  ##  do this again, since the order may have changed
  num.DifExpMeasures <- c(rep(0,length(DifExp.list)))
  for(i in 1:num.studies)
  {
    temp.label <- (DifExp.list[[i]]$label)
    t <- as.factor(levels(as.factor(temp.label)))
    num.DifExpMeasures[i] <- length(t)-1
  }

  ##  check the newnames data frame for bad entries
  starting.newnames.length <- length(newnames[,1])
  temp.newnames <- newnames[!is.na(newnames$new.name),]
  newnames <- temp.newnames
  temp.newnames <- newnames[!is.na(newnames$old.name),]
  newnames <- temp.newnames
  temp.newnames <- newnames[!is.na(newnames$chip),]
  newnames <- temp.newnames
  final.newnames.length <- length(newnames[,1])
  if (final.newnames.length != starting.newnames.length)
  {
    if (show.warnings)
      cat(starting.newnames.length - final.newnames.length, "rows removed from newnames data frame due to NA values\n")
  }

  newnames.char <- cbind(as.character(newnames[,1]),as.character(newnames[,2]),as.character(newnames[,3]))
  chipset.names.char <- newnames.char[,1]
  old.names.char <- newnames.char[,2]
  new.names.char <- newnames.char[,3] 
  num.genes <- as.integer(length(newnames.char[,1]))

  names.only.list <- list()
  for(i in 1:num.studies)
  {  names.only.list[[i]] <- as.vector(DifExp.list[[i]][,2])  }

  study.gene.names <- list()
  for(i in 1:num.studies)
  {  study.gene.names[[i]] <- as.vector(DifExp.list[[i]][,2])  }

  trimmed.difexp.list <- list()
  for(i in 1:num.studies)
  {
    temp.study <- DifExp.list[[i]]
    temp.width <- length(temp.study[1,])
    trimmed.difexp.list[[i]] <- cbind(as.matrix(temp.study[3:temp.width]))
  }

  num.genes.per.study <- rep(0, num.studies)
  for(i in 1:num.studies)
    num.genes.per.study[i] <- length(DifExp.list[[i]][,1])
  num.genes.per.study <- as.integer(num.genes.per.study)

  cov.mat.offsets <- rep(0, num.studies)
  for (i in 2:num.studies)
    cov.mat.offsets[i] <- sum(num.DifExpMeasures[1:(i-1)])

  ##study.name.sort.index <- .Call("sort_name_lists", names.only.list, PACKAGE="metahdep.dll")
  ##newnames.sort.index <- .Call("sort_newname_list", new.names.char, PACKAGE="metahdep.dll")
  study.name.sort.index <- .Call("sort_name_lists", names.only.list)
  newnames.sort.index <- .Call("sort_newname_list", new.names.char)

  return.list <- list()

  gene.list <- as.vector(unique(newnames$new.name))
  gene.list.length <- length(gene.list)
  update.step <- trunc(gene.list.length/5)

  current.list.index <- 1

  metah.update.indices <- trunc(seq(1, gene.list.length, length.out=22))
  cat("\n|-----Formatting-----|\n")

  ###########################################################
  ##	BEGIN CYCLING THROUGH THE GENE LIST
  ###########################################################

  for (gene_id in 1:gene.list.length)
  {
    if (is.element(gene_id, metah.update.indices))
    {
      cat("*")
      flush.console()
    }

    current.gene <- gene.list[gene_id]

    ##old.name.indices2 <- .Call("find_old_names2", new.names.char, newnames.sort.index, current.gene, PACKAGE="metahdep.dll")
    old.name.indices2 <- .Call("find_old_names2", new.names.char, newnames.sort.index, current.gene)
    ##  only need to continue processing the gene if old.names were found
    if (is.null(old.name.indices2))
    {
      if (show.warnings)
        cat("No old.name matches found for newname: ", current.gene, "\n")
    }
    else
    {
      old.name.list <- old.names.char[old.name.indices2+1]
      old.chipsets <- chipset.names.char[old.name.indices2+1]

      ##row.indices <- .Call("get_row_indices2", study.gene.names, old.name.list, old.chipsets, full.chipset.list, study.name.sort.index, PACKAGE="metahdep.dll")
      row.indices <- .Call("get_row_indices2", study.gene.names, old.name.list, old.chipsets, full.chipset.list, study.name.sort.index)
      ##  only need to continue processing the gene if row.indices were found
      if (is.null(row.indices))
      {
        if (show.warnings == TRUE)
          cat("No data found (in any study) for gene: ", current.gene, "(", gene_id, ")\n")
      }
      else
      {
        colnames(row.indices) <- c("Study", "Row")
        row.indices <- row.indices+1

        current.study.list <- row.indices[,1]
        current.row.list <- row.indices[,2]
        total.difexpmeasures <- sum(num.DifExpMeasures[current.study.list])
        num.rows.found <- length(current.row.list)

        theta <- rep(0.0, total.difexpmeasures)
        V <- matrix(0.0, total.difexpmeasures, total.difexpmeasures)
        X <- matrix(0.0, total.difexpmeasures, p)
        M <- matrix(0.0, total.difexpmeasures, total.difexpmeasures)
        block.offset <- 0

        #####################################################
        ##   BUILD THETA AND V
        #####################################################
        for (i in 1:num.rows.found)
        {
          current.study <- current.study.list[i]
          current.row <- current.row.list[i]
          current.num.difexpmeasures <- num.DifExpMeasures[current.study]
          current.data <- trimmed.difexp.list[[current.study]][current.row,]

          ##  Append the DifExpMeasures to theta
          for (j in 1:current.num.difexpmeasures)
          {  theta[block.offset + j] = current.data[j]  }

          ##  Add the new var/cov values to V  (as a block on the diagonal)
          varcov.offset <- current.num.difexpmeasures
          for (j in 1:current.num.difexpmeasures)
          {
            for (k in j:current.num.difexpmeasures)
            {
              V[block.offset + j, block.offset + k] = current.data[varcov.offset + k]
              if (j == k)
              {
                if (current.data[varcov.offset + k] < min.var)
                {  V[block.offset + j, block.offset + k] = min.var  }
              }
              else
              {  V[block.offset + k, block.offset + j] = current.data[varcov.offset + k]  }
            }
            varcov.offset <- varcov.offset + current.num.difexpmeasures - j
          }
          block.offset <- block.offset + current.num.difexpmeasures
        }

        ##  if pd.verify == TRUE then check V to verify that it is positive definite
        ##  this should only happen if one of the data frames in the DifExp.list is formatted incorrectly
        ##  or if one of them simply has bad data for the variance/covariances
        if (pd.verify)
        {
          ##  do a cholesky decomposition.  this will fail if V is not positive definite
          ##  on success, pd.check is a matrix.  on failure, pd.check is not a matrix
          pd.check <- try(chol(V))
          if (!is.matrix(pd.check))
          {
            cat("The variance/covariance matrix for gene", current.gene, "(",gene_id,")", "is not positive definite\n")
            cat("Study / Row indices:\n")
            for (i in 1:nrow(row.indices))
              cat(row.indices[i,], "\n")
            flush.console()
          }
        }

        cov.mat.row.indices <- rep(0, total.difexpmeasures)
        cov.mat.row.indices.offset <- 0
        #####################################################
        ##   BUILD X
        #####################################################
        ##  check for the case of a 1x1 design matrix
        if (total.difexpmeasures == 1)
        {
          X <- as.matrix(1)
          colnames(X) <- as.vector(covariate.names[1])
        }
        else
        {
          for (i in 1:num.rows.found)
          {
            for (j in 1:num.DifExpMeasures[current.study.list[i]])
              cov.mat.row.indices[cov.mat.row.indices.offset + j] <- cov.mat.offsets[current.study.list[i]] + j
            cov.mat.row.indices.offset <- cov.mat.row.indices.offset + num.DifExpMeasures[current.study.list[i]]
          }
          X <- covariates[cov.mat.row.indices,]

          ##   CHECK THAT X IS FULL-RANK
          ##   (code copied from original JRS implementation)
          X.0 <- matrix(X[,1],ncol=1)
          X.names <- covariate.names[1]
          for(i in 2:p)
          {
            X.temp <- cbind(X.0,X[,i])
            if(qr(X.temp)$rank==ncol(X.temp))
            {
              X.0 <- X.temp
              X.names <- c(X.names,covariate.names[i])
            }
          }
          X.revised <- X.0
          colnames(X.revised) <- X.names
          X <- X.revised
        }

        #####################################################
        ##   BUILD M
        #####################################################
        max.k <- 1
        block.offset <- 0

        ##  if all observed differential expression measures come from the same dependence group, or
        ##  if all observations come from different dependence groups (no 2 elements in theta are in the
        ##  same dependence group) then it doesn't make sense to perform delta-splitting for this gene.
        ##  In this case, M should be all 0's and max.k should also be 0
        current.num.dep.groups <- length(unique(dependence.groups[current.study.list]))
        if ( current.num.dep.groups == 1 | current.num.dep.groups == total.difexpmeasures)
        {
          M <- matrix(0, total.difexpmeasures, total.difexpmeasures)
          max.k <- 0
        }
        else
        {
          for (i in 1:num.dep.groups)
          {
            ##  find how big the dependence group block needs to be
            num.dep.rows <- 0
            for (j in 1:num.rows.found)
            {
              if (dependence.groups[current.study.list[j]] == i)
                num.dep.rows <- num.dep.rows + num.DifExpMeasures[current.study.list[j]]
            }
            if (num.dep.rows > max.k)
              max.k <- num.dep.rows

            ##  add the block for the current dependence group to the diagonal of M, if needed
            if (num.dep.rows > 0)
            {
              for (j in 1:num.dep.rows)
              {
                for (k in 1:num.dep.rows)
                {
                  M[block.offset + j, block.offset + k] <- 1
                  M[block.offset + k, block.offset + j] <- 1
                }
              }
            }
            block.offset <- block.offset + num.dep.rows
          }
          ##  re-set the diagonal elements to 0
          for (i in 1:total.difexpmeasures)
            M[i,i] <- 0
        }

        ##  if V checks out, then add everything to the prep list
        if (is.matrix(pd.check))
        {
          ##  add theta, V, X, and M to the return list.  add the gene name on the end for reference
          ##  if include.row.indices is TRUE, then return those also
          if (include.row.indices == TRUE)
            return.list[[current.list.index]] <- new("metaprep", theta=theta, V=V, X=X, M=M, max.k=as.integer(max.k), row.indices=row.indices, gene=current.gene)
          else
            return.list[[current.list.index]] <- new("metaprep", theta=theta, V=V, X=X, M=M, max.k=as.integer(max.k), gene=current.gene)

          ##  increment the index into the return.list
          current.list.index <- current.list.index + 1
        }
      }  ##  if (!is.null(row.indices))
    }  ##  if (!is.null(old.names)
  }  ##  for (gene_id in 1:gene.list.length)

  num.total.genes.found <- length(return.list)

  cat("\nCreated entries for", num.total.genes.found, "genes\n")
  return(return.list)
}


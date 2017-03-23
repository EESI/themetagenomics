jsd <- function(p,q) {

  m <- .5*(p + q)
  div <- .5*sum(p*log(p/m)) + .5*sum(q*log(q/m))
  dis <- sqrt(div)

  return(dis)
}

norm10 <- function(x) (x-min(x))/(max(x)-min(x))


create_modelframe <- function(formula,refs,metadata){

  rnames <- rownames(metadata)

  metadata <- as.data.frame(unclass(model.frame(formula,data=metadata)))
  term_names <- colnames(metadata)

  cnames <- NULL
  modelframe <- NULL
  for (i in seq_along(term_names)){

    term <- term_names[i]

    if (class(metadata[,term]) != 'factor'){
      tmp <- metadata[term]
      tmp <- model.matrix(as.formula(paste0('~',term)),tmp)
      cnames <- c(cnames,colnames(tmp)[-1])
      modelframe <- cbind(modelframe,tmp[,-1])
    }else{
      tmp <- metadata[term]
      tmp[[1]] <- relevel(tmp[[1]],ref=refs[1])
      refs <- refs[-1]
      mm <- model.matrix(as.formula(paste0('~',term)),tmp)
      cnames <- c(cnames,colnames(mm)[-1])

      modelframe <- data.frame(modelframe,
                               lapply(seq_along(levels(tmp[[1]]))[-1],
                                      function(j){
                                        relevel(as.factor(ifelse(mm[,j] == 1,levels(tmp[[1]])[j],levels(tmp[[1]])[1])),ref=levels(tmp[[1]])[1])
                                      }
                               )
      )

    }

  }

  colnames(modelframe) <- cnames
  rownames(modelframe) <- rnames

  return(modelframe)

}


make_ppd_x <- function(covariate,mod,modelframe_full,npoints){

  if (class(modelframe_full[[covariate]]) == 'numeric'){

    if (missing(mod)){
      x <- matrix(0.0,nrow=npoints,ncol=ncol(modelframe_full) + 1,dimnames=list(NULL,c('Intercept',colnames(modelframe_full))))
      x[,1] <- 1
      for (i in 1:ncol(modelframe_full)){
        if (colnames(modelframe_full)[i] == covariate){
          x[,i+1] <- seq(min(modelframe_full[[i]]),max(modelframe_full[[i]]),length.out=npoints)
        }else if (colnames(modelframe_full)[i] == 'numeric'){
          x[,i+1] <- median(modelframe_full[[i]])
        }else{
          x[,i+1] <- sum(levels(modelframe_full[[i]])[2] == modelframe_full[[i]])/length(modelframe_full[[i]])
        }
      }

      out <- list(x)
    }else{
      x <- matrix(0.0,nrow=npoints,ncol=ncol(modelframe_full) + 1,dimnames=list(NULL,c('Intercept',colnames(modelframe_full))))
      x[,1] <- 1
      for (i in 1:ncol(modelframe_full)){
        if (colnames(modelframe_full)[i] == covariate){
          x[,i+1] <- seq(min(modelframe_full[[i]]),max(modelframe_full[[i]]),length.out=npoints)
        }else if (colnames(modelframe_full)[i] == mod){
          x[,i+1] <- 1
        }else if (colnames(modelframe_full)[i] == 'numeric'){
          x[,i+1] <- median(modelframe_full[[i]])
        }else{
          x[,i+1] <- sum(levels(modelframe_full[[i]])[2] == modelframe_full[[i]])/length(modelframe_full[[i]])
        }
      }
      x1 <- x
      x2 <- x
      x2[,mod] <- 0

      out <- list(x1,x2)
      names(out) <- c(mod,paste0('not',mod))
      out <- list(out)
    }

    names(out) <- covariate

    return(out)

  }

  if (class(modelframe_full[[covariate]]) == 'factor'){

    x <- matrix(0.0,nrow=2,ncol=ncol(modelframe_full) + 1,dimnames=list(NULL,c('Intercept',colnames(modelframe_full))))
    x[,1] <- 1
    for (i in 1:ncol(modelframe_full)){
      if (colnames(modelframe_full)[i] == covariate){
        x[,i+1] <- c(0,1)
      }else if (colnames(modelframe_full)[i] == 'numeric'){
        x[,i+1] <- median(modelframe_full[[i]])
      }else{
        x[,i+1] <- sum(levels(modelframe_full[[i]])[2] == modelframe_full[[i]])/length(modelframe_full[[i]])
      }
    }

    out <- list(x)
    names(out) <- covariate

    return(out)

  }

}

create_multiclasses_table <- function(modelframe,modelframe_full){
  classes <- sapply(modelframe,class)
  multiclasses <- classes
  multiclasses[which(sapply(modelframe,function(x) length(levels(x))) > 2)] <- 'multiclass'
  multiclasses <- sapply(seq_along(multiclasses), function(i) {

    if (classes[i] == 'factor') {
      j <- length(levels(modelframe[,i]))-1
      class <- rep(classes[i],j)
      multiclass <- rep(multiclasses[i],j)
      lev <- levels(modelframe[,i])[-1]
      ref <- levels(modelframe[,i])[1]
    }else{
      j <- 1
      class <- rep(classes[i],j)
      multiclass <- rep(multiclasses[i],j)
      lev <- 2
      ref <- 1
    }

    cbind(class,multiclass,lev,ref)

  })
  cbind(full=names(sapply(modelframe_full,class)),
        do.call('rbind',multiclasses))
}

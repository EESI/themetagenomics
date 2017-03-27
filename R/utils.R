jsd <- function(p,q) {

  m <- .5*(p + q)
  div <- .5*sum(p*log(p/m)) + .5*sum(q*log(q/m))
  dis <- sqrt(div)

  return(dis)
}

norm10 <- function(x) (x-min(x))/(max(x)-min(x))

normz <- function(x) (x-mean(x))/sd(x)

# dplyr dense_rank
dense_rank <- function(x) {

  r <- rank(x,na.last='keep')
  match(r,sort(unique(r)))

}

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


make_ppd_x <- function(estimated_effects,covariate,mod,npoints=100){

  formula <- estimated_effects$formula
  metadata <- estimated_effects$data
  modelframe_full <- estimated_effects$modelframe_full

  if (missing(covariate)){
    if (check_for_splines(formula,metadata)){
      modelframe <- model.matrix(formula,data=metadata)
      ppd_idx <- which(labels(modelframe)[[2]] %in% colnames(modelframe_full))
      modelframe[,ppd_idx] <- matrix(rep(colMeans(modelframe[,ppd_idx]),nrow(modelframe)),nrow(modelframe),byrow=TRUE)
      return(list(modelframe))
    }else{
      modelframe <- model.matrix(formula,data=metadata)
      return(list(modelframe))
    }
  }

  if (attr(covariate,'baseclass') == 'factor'){
    modelframe <- colMeans(model.matrix(formula,data=metadata))
    modelframe <- rbind(modelframe,modelframe)
    modelframe[1,covariate] <- 0
    modelframe[2,covariate] <- 1
    return(list(modelframe))
  }

  if (attr(covariate,'multiclass') == 'spline'){
    newdata <- metadata
    newdata[[covariate]] <- seq(min(newdata[[covariate]]),max(newdata[[covariate]]),length.out=nrow(newdata))
    modelframe <- model.matrix(formula,data=newdata)
    ppd_idx <- which(labels(modelframe)[[2]] %in% colnames(modelframe_full))
    modelframe[,ppd_idx] <- matrix(rep(colMeans(modelframe[,ppd_idx]),nrow(modelframe)),nrow(modelframe),byrow=TRUE)
  }

  if (attr(covariate,'baseclass') == 'numeric' & attr(covariate,'multiclass') != 'spline'){
    newdata <- metadata
    newdata[[covariate]] <- seq(min(newdata[[covariate]]),max(newdata[[covariate]]),length.out=nrow(newdata))
    modelframe <- model.matrix(formula,data=newdata)
    modelframe[,colnames(modelframe) != covariate] <- matrix(rep(colMeans(modelframe[,colnames(modelframe) != covariate]),
                                                                 nrow(modelframe)),nrow(modelframe),byrow=TRUE)
  }

  if (!missing(mod)){
    mf0 <- mf1 <- modelframe
    mf0[,mod] <- 0
    mf1[,mod] <- 1
    modelframe <- list(mf0=mf0,mf1=mf1)
  }else{
    modelframe <- list(modelframe)
  }
  return(modelframe)

}

create_multiclasses_table <- function(modelframe,modelframe_full,splines=NULL){
  classes <- sapply(modelframe,class)
  multiclasses <- classes
  multiclasses[which(sapply(modelframe,function(x) length(levels(x))) > 2)] <- 'multiclass'
  multiclasses[splines] <- 'spline'
  multiclasses <- sapply(seq_along(multiclasses), function(i) {

    if (classes[i] == 'factor') {
      j <- length(levels(modelframe[,i]))-1
      baseclass <- rep(classes[i],j)
      multiclass <- rep(multiclasses[i],j)
      lev <- levels(modelframe[,i])[-1]
      ref <- levels(modelframe[,i])[1]
    }else{
      j <- 1
      baseclass <- rep(classes[i],j)
      multiclass <- rep(multiclasses[i],j)
      lev <- 2
      ref <- 1
    }

    cbind(baseclass,multiclass,lev,ref)

  })

  out <- cbind(full=names(sapply(modelframe_full,class)),
        do.call('rbind',multiclasses))
  rownames(out) <- NULL
  out <- data.frame(out,stringsAsFactors=FALSE)

  return(out)

}

check_for_splines <- function(formula,metadata){

  vars <- terms(formula,data=metadata,
                specials=c('s','bs','ns','poly'))

  if (any(!sapply(attr(vars,'specials'),is.null))) return(TRUE) else return(FALSE)

}

extract_spline_info <- function(formula,metadata,remove_only=FALSE){

  vars <- terms(formula,data=metadata,
                specials=c('s','bs','ns','poly'))
  splines <- attr(vars,'special')
  splines <- splines[!sapply(splines,is.null)]

  vars_old <- labels(vars)
  vars_new <- all.vars(formula,unique=FALSE)
  rhs <- as.character(formula)[2]

  for (i in seq_along(vars_old)){
    rhs <- gsub(vars_old[i],vars_new[i],rhs,fixed=TRUE)
  }

  formula_new <- as.formula(paste0('~',rhs))

  if (remove_only) return(formula_new)

  info <- vector(mode='list',length=length(unlist(splines)))
  j <- 1
  for (s in names(splines)){
    for (i in seq_along(splines[[s]])){
      info[[j]] <- list(
        var=vars_new[i],
        spline=s,
        expansion=model.frame(as.formula(paste0('~',vars_old[i])),data=metadata)
      )
      j <- j + 1
    }
  }

  return(list(info=info,formula=formula_new))

}

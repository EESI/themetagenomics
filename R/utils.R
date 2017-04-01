# Jensen Shannon distance
jsd <- function(p,q) {

  m <- .5*(p + q)
  div <- .5*sum(p*log(p/m)) + .5*sum(q*log(q/m))
  dis <- sqrt(div)

  return(dis)
}

# Normalize to 0-1 range
norm10 <- function(x) (x-min(x))/(max(x)-min(x))

# Z score normalize
normz <- function(x) (x-mean(x))/stats::sd(x)

# dplyr dense_rank
dense_rank <- function(x) {

  r <- rank(x,na.last='keep')
  match(r,sort(unique(r)))

}

# Create modelframe that expands catagoricals into dummy variables
create_modelframe <- function(formula,metadata,refs){

  rnames <- rownames(metadata)

  metadata <- as.data.frame(unclass(stats::model.frame(formula,data=metadata)))
  term_names <- colnames(metadata)

  cnames <- NULL
  modelframe <- list()
  for (i in seq_along(term_names)){

    term <- term_names[i]

    if (class(metadata[,term]) != 'factor'){
      tmp <- metadata[term]
      tmp <- stats::model.matrix(stats::as.formula(paste0('~',term)),tmp)
      cnames <- c(cnames,colnames(tmp)[-1])
      modelframe <- c(modelframe,list(tmp[,-1]))
    }else{
      tmp <- metadata[term]
      tmp[[1]] <- stats::relevel(tmp[[1]],ref=refs[1])
      refs <- refs[-1]
      mm <- stats::model.matrix(stats::as.formula(paste0('~',term)),tmp)
      cnames <- c(cnames,colnames(mm)[-1])

      modelframe <- c(modelframe,lapply(seq_along(levels(tmp[[1]]))[-1],
                           function(j){
                             stats::relevel(as.factor(ifelse(mm[,j] == 1,levels(tmp[[1]])[j],levels(tmp[[1]])[1])),
                                            ref=levels(tmp[[1]])[1])
                           }))

    }

  }

  modelframe <- as.data.frame(modelframe)
  colnames(modelframe) <- cnames
  rownames(modelframe) <- rnames

  return(modelframe)

}

# Create design matrix for ppd sampling
make_ppd_x <- function(estimated_effects,covariate,mod,npoints=100){

  formula <- estimated_effects$formula
  metadata <- estimated_effects$data
  modelframe_full <- estimated_effects$modelframe_full

  if (missing(covariate)){
    if (check_for_splines(formula,metadata)){
      modelframe <- stats::model.matrix(formula,data=metadata)
      ppd_idx <- which(labels(modelframe)[[2]] %in% colnames(modelframe_full))
      modelframe[,ppd_idx] <- matrix(rep(colMeans(modelframe[,ppd_idx]),nrow(modelframe)),nrow(modelframe),byrow=TRUE)
      return(list(modelframe))
    }else{
      modelframe <- stats::model.matrix(formula,data=metadata)
      return(list(modelframe))
    }
  }

  if (attr(covariate,'baseclass') == 'factor'){
    modelframe <- colMeans(stats::model.matrix(formula,data=metadata))
    modelframe <- rbind(modelframe,modelframe)
    modelframe[1,covariate] <- 0
    modelframe[2,covariate] <- 1
    return(list(modelframe))
  }

  if (attr(covariate,'multiclass') == 'spline'){
    newdata <- metadata
    newdata[[covariate]] <- seq(min(newdata[[covariate]]),max(newdata[[covariate]]),length.out=nrow(newdata))
    modelframe <- stats::model.matrix(formula,data=newdata)
    ppd_idx <- which(labels(modelframe)[[2]] %in% colnames(modelframe_full))
    modelframe[,ppd_idx] <- matrix(rep(colMeans(modelframe[,ppd_idx]),nrow(modelframe)),nrow(modelframe),byrow=TRUE)
  }

  if (attr(covariate,'baseclass') == 'numeric' & attr(covariate,'multiclass') != 'spline'){
    newdata <- metadata
    newdata[[covariate]] <- seq(min(newdata[[covariate]]),max(newdata[[covariate]]),length.out=nrow(newdata))
    modelframe <- stats::model.matrix(formula,data=newdata)
    if (ncol(modelframe) > 2)
      modelframe[,colnames(modelframe) != covariate,drop=FALSE] <- matrix(rep(colMeans(modelframe[,colnames(modelframe) != covariate,drop=FALSE]),
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

# Create a lookup table that identifies splines and multiclass factors
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

  if (class(multiclasses) == 'list'){
    out <- cbind(full=names(sapply(modelframe_full,class)),
                 do.call('rbind',multiclasses))
    rownames(out) <- NULL
  }else{
    out <- cbind(colnames(modelframe_full),t(multiclasses))
    colnames(out) <- c('full','baseclass','multiclass','lev','ref')
    rownames(out) <- NULL
  }

  out <- data.frame(out,stringsAsFactors=FALSE)

  return(out)

}

# Checks whether a continuous covariate involves a basis function
check_for_splines <- function(formula,metadata){

  vars <- stats::terms(formula,data=metadata,
                specials=c('s','bs','ns','poly'))

  if (any(!sapply(attr(vars,'specials'),is.null))) return(TRUE) else return(FALSE)

}

# Extract basis function design matrix and basis type
extract_spline_info <- function(formula,metadata,remove_only=FALSE){

  vars <- stats::terms(formula,data=metadata,
                specials=c('s','bs','ns','poly'))
  splines <- attr(vars,'special')
  splines <- splines[!sapply(splines,is.null)]

  vars_old <- labels(vars)
  vars_new <- all.vars(formula,unique=FALSE)
  rhs <- as.character(formula)[2]

  for (i in seq_along(vars_old)){
    rhs <- gsub(vars_old[i],vars_new[i],rhs,fixed=TRUE)
  }

  formula_new <- stats::as.formula(paste0('~',rhs))

  if (remove_only) return(formula_new)

  info <- vector(mode='list',length=length(unlist(splines)))
  j <- 1
  for (b in names(splines)){
    for (i in seq_along(splines[[b]])){
      info[[j]] <- list(var=vars_new[splines[[b]][i]],
                        spline=b,
                        expansion=stats::model.frame(stats::as.formula(paste0('~',vars_old[splines[[b]][i]])),data=metadata))
      j <- j + 1
    }
  }

  return(list(info=info,formula=formula_new))

}

# create_skeleton from rstan
create_skel <- function (pars,dims){
  lst <- lapply(seq_along(pars),function(i){
    len_dims <- length(dims[[i]])
    if (len_dims < 1) return(0)
    return(array(0, dim = dims[[i]]))
  })
  names(lst) <- pars
  lst
}

# rmvnorm from stm
rmvnorm <- function(n,mu,Sigma,chol_Sigma=chol(Sigma)){
  E <- matrix(rnorm(n*length(mu)),n,length(mu))
  t(t(E %*% chol_Sigma) + c(mu))
}

# simBetas from stm
ppd_weights <- function(parameters,nsims=100){
  betas <- vector(mode='list',length=length(parameters))
  for (i in seq_along(parameters)){
    betas[[i]] <- do.call('rbind',
                          lapply(parameters[[i]],function(x) rmvnorm(n=nsims,mu=x$est,Sigma=x$vcov)))
  }
  return(betas)
}

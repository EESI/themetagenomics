#' Estimate topic effects for a given covariate
#'
#' Given a covariate of interest, measure its relationship with the samples over topics distribution estimated by the STM.
#'
#' @param covariate String corresponding to the name of the covariate of interest found in metadata.
#' @param topics STM object containing the topic information
#' @param metadata Dataframe or matrix of metadata. If the STM was fit with covariate information, this must point to the same metadata used in the STM.
#' @param nsims (optional) Number of simulations to perform for inferring covariate effects.
#' @param ui_level With of uncertainty interval to report effects. Defaults to .95.
#' @return A list containing
#'
#' \item{est}{The estimate of the effect}
#' \item{rank}{The rank of the topics based on the effect estimate (increasing)}
#' \item{sig}{A vector of topic indexes for topics whose uncertainty intervals did not enclose 0.}
#' @export

estimate_topic_effects <- function(topics,metadata,formula,refs,moderator,nsims=100,ui_level=.8,...){

  fit <- topics$fit
  K <- fit$settings$dim$K

  if (missing(formula)){
    formula <- fit$settings$covariates$formula
    if (is.null(formula)) stop('Please provide a formula.')
  }

  modelframe <- as.data.frame(unclass(model.frame(formula,metadata)))
  classes <- sapply(modelframe,class)
  multiclasses <- classes
  multiclasses[which(sapply(modelframe,function(x) length(levels(x))) > 2)] <- 'multiclass'

  if (sum(classes == 'factor') > 0){
    if (missing(refs)){
      warning('References are recommended for factors. Using the first level(s).')
      refs <- unlist(lapply(modelframe[,classes == 'factor'],function(x) levels(x)[1]))
    }

    if (sum(classes == 'factor') != length(refs)) stop('A reference is required for each factor.')
    ref_check <- all(sapply(seq_along(refs), function(i) refs[i] %in% lapply(modelframe[,classes == 'factor'],levels)[[i]]))
    if (!ref_check) stop('Reference(s) not found in factor(s).')

    j <- 1
    for (i in seq_along(classes)){
      if (classes[i] == 'factor'){
        modelframe[,i] <- relevel(as.factor(modelframe[,i]),ref=refs[j])
        j <- j+1
      }
    }

  }

  if (missing(moderator)){
    if (sum(multiclasses == 'factor') == 1){
      moderator <- colnames(modelframe)[multiclasses == 'factor']
    }else if (sum(multiclasses == 'factor') > 1){
      moderator <- colnames(modelframe)[multiclasses == 'factor'][1]
      warning(sprintf('No moderator specified. Setting it to %s.',moderator))
    }else{
      moderator <- NULL
    }
  }

  formula <- as.formula(sprintf('1:%s %s',K,paste0(formula,collapse=' ')))
  estimated_effects <- stm::estimateEffect(formula,fit,modelframe,uncertainty='Global')

  estimated_effects$modelframe_full <- topics$modelframe
  estimated_effects$modelframe <- modelframe
  estimated_effects$moderator <- moderator

  topic_effects <- estimate_topic_effects_backend(estimated_effects,nsims=nsims,ui_level=ui_level,...)

  return(list(topic_effects=topic_effects,modelframe=topics$modelframe))

}

#  Backend to extract effects for \code{\link{estimate_topic_effects}}.

estimate_topic_effects_backend <- function(estimated_effects,nsims=100,ui_level=.8,npoints=100){

  K <- length(estimated_effects$topics)

  ui_offset <- (1-ui_level)/2
  ui_interval <- paste0(c(ui_offset,1-ui_offset)*100,'%')
  ui_levels <- matrix(0.0,K,2,dimnames=list(NULL,ui_interval))

  modelframe <- estimated_effects$modelframe

  cov_switch <- estimated_effects$moderator
  cov_switch_levels <- levels(modelframe[,cov_switch])

  simbetas <- stm:::simBetas(estimated_effects$parameters,nsims=nsims)
  for (i in seq_along(simbetas)) colnames(simbetas[[i]]) <- names(estimated_effects$parameter[[1]][[1]]$est)

  # covariate_list <- vector(mode='list',length=ncol(simbetas[[1]])-1)
  # names(covariate_list) <- colnames(simbetas[[1]])[-1]

  covariate_list <- vector(mode='list',length=ncol(estimated_effects$modelframe_full))
  names(covariate_list) <- colnames(estimated_effects$modelframe_full)

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
  multiclasses <- do.call('rbind',multiclasses)

  for (i in seq_len(nrow(multiclasses))){

    # freezes other covariates except target

    cov_i <- rownames(multiclasses)[i]

    if (multiclasses[i,'class'] == 'numeric'){

      cthis <- stm:::produce_cmatrix(prep=estimated_effects,covariate=cov_i,
                                     method='continuous',npoints=npoints,
                                     moderator=cov_switch,moderator.value=cov_switch_levels[2])

      cthis_switch <- stm:::produce_cmatrix(prep=estimated_effects,covariate=cov_i,
                                     method='continuous',npoints=npoints,
                                     moderator=cov_switch,moderator.value=cov_switch_levels[1])

      cdata_switch <- cthis_switch$cdata
      cmat_switch <- cthis_switch$cmatrix

    }else if (multiclasses[i,'class'] == 'factor'){

      factor_levels <- levels(modelframe[,cov_i])
      cthis <- stm:::produce_cmatrix(prep=estimated_effects,covariate=cov_i,
                                     method='difference',
                                     cov.value1=multiclasses[i,'ref'],cov.value2=multiclasses[i,'lev'])

    }

    cdata <- cthis$cdata
    cmat <- cthis$cmatrix

    fitted <- NA
    fitted_switch <- NA
    if (multiclasses[i,'class'] == 'numeric'){

      fitted <- lapply(seq_len(K),function(x) matrix(0.0,npoints,4,dimnames=list(NULL,c('estimate',ui_interval,'covariate'))))
      names(fitted) <- paste0('T',1:K)
      fitted_switch <- lapply(seq_len(K),function(x) matrix(0.0,npoints,4,dimnames=list(NULL,c('estimate',ui_interval,'covariate'))))
      names(fitted_switch) <- paste0('T',1:K)

    }

    est <- lapply(seq_len(K), function(x) matrix(0.0,1,3,dimnames=list('slope',c('estimate',ui_interval))))
    for (j in seq_len(K)){

      sims <- cmat %*% t(simbetas[[j]])

      if (multiclasses[i,'class'] == 'numeric'){

        sims_switch <- cmat_switch %*% t(simbetas[[j]])

        est[[j]][,] <- c(mean(simbetas[[j]][,cov_i]),quantile(simbetas[[j]][,cov_i],c(ui_offset,1-ui_offset)))

        fitted[[j]][,] <- c(rowMeans(sims),t(apply(sims,1,function(x) quantile(x,c(ui_offset,1-ui_offset)))),cdata[,cov_i])
        fitted_switch[[j]][,] <- c(rowMeans(sims_switch),t(apply(sims_switch,1,function(x) quantile(x,c(ui_offset,1-ui_offset)))),cdata_switch[,cov_i])

      }else if (multiclasses[i,'class'] == 'factor'){

        diff <- sims[2,] - sims[1,]
        est[[j]][,] <- c(mean(diff),quantile(diff,c(ui_offset,1-ui_offset)))

      }

    }

    est_mat <- do.call('rbind',est)
    rownames(est_mat) <- paste0('T',1:K)
    rank <- dense_rank(est_mat[,1])
    names(rank) <- paste0('T',1:K)
    sig <- which(rowSums(sign(est_mat[,2:3])) != 0)

    covariate_list[[i]] <- list(est=est_mat,rank=rank,sig=sig,fitted=fitted,fitted_switch=fitted_switch)

  }

  return(covariate_list)

}

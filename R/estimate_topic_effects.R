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

estimate_topic_effects <- function(topics,metadata,formula,nsims=100,ui_level=.8,...){

  fit <- topics$fit
  K <- fit$settings$dim$K

  if (missing(formula)){
    formula <- fit$settings$covariates$formula
    if (is.null(formula)) stop('Please provide a formula.')
  }

  formula <- as.formula(sprintf('1:%s %s',K,paste0(formula,collapse=' ')))

  estimated_effects <- stm::estimateEffect(formula,fit,metadata,uncertainty='Global')

  modelframe <- as.data.frame(unclass(estimated_effects$modelframe))
  estimated_effects$modelframe <- modelframe

  # must provide binary convariate as either factor or character, even if 0-1, otherwise treated as continuous

  modelframe_classes <- sapply(modelframe,class)
  if (NCOL(modelframe) == 0){
    stop('modelframe has 0 dimensions.')
  }else if (any(sapply(modelframe,function(x) length(levels(x))) > 2)){
    stop('Recode multilevel factor as dummy variables.')
  }else{
    topic_effects <- estimate_topic_effects_backend(estimated_effects,nsims=nsims,ui_level=ui_level,...)
  }

  return(topic_effects)

}

#  Backend to extract effects for \code{\link{estimate_topic_effects}}.

estimate_topic_effects_backend <- function(estimated_effects,nsims=100,ui_level=.8,npoints=100){

  K <- length(estimated_effects$topics)

  ui_offset <- (1-ui_level)/2
  ui_interval <- paste0(c(ui_offset,1-ui_offset)*100,'%')
  ui_levels <- matrix(0.0,K,2,dimnames=list(NULL,ui_interval))

  modelframe <- estimated_effects$modelframe
  modelframe_classes <- sapply(modelframe,class)
  cov_switch <- names(which(modelframe_classes == 'factor'))
  cov_switch_levels <- levels(modelframe[,cov_switch])

  simbetas <- stm:::simBetas(estimated_effects$parameters,nsims=nsims)
  for (i in seq_along(simbetas)) colnames(simbetas[[i]]) <- c('Intercept',names(modelframe_classes))


  covariate_list <- vector(mode='list',length=length(modelframe_classes))
  names(covariate_list) <- names(modelframe_classes)
  for (i in seq_along(modelframe_classes)){

    # freezes other covariates except target

    cov_i <- names(modelframe_classes)[i]

    if (modelframe_classes[i] == 'numeric'){

      cthis <- stm:::produce_cmatrix(prep=estimated_effects,covariate=cov_i,
                                     method='continuous',npoints=npoints,
                                     moderator=cov_switch,moderator.value=cov_switch_levels[2])

      cthis_switch <- stm:::produce_cmatrix(prep=estimated_effects,covariate=cov_i,
                                     method='continuous',npoints=npoints,
                                     moderator=cov_switch,moderator.value=cov_switch_levels[1])

      cdata_switch <- cthis_switch$cdata
      cmat_switch <- cthis_switch$cmatrix

    }else if (modelframe_classes[i] == 'factor'){

      factor_levels <- levels(modelframe[,cov_i])
      cthis <- stm:::produce_cmatrix(prep=estimated_effects,covariate=cov_i,
                                     method='difference',
                                     cov.value1=factor_levels[1],cov.value2=factor_levels[2])

    }

    cdata <- cthis$cdata
    cmat <- cthis$cmatrix

    fitted <- NA
    fitted_switch <- NA
    if (modelframe_classes[i] == 'numeric'){
      fitted <- lapply(seq_len(K),function(x) matrix(0.0,npoints,4,dimnames=list(NULL,c('estimate',ui_interval,'covariate'))))
      names(fitted) <- paste0('T',1:K)
      fitted_switch <- lapply(seq_len(K),function(x) matrix(0.0,npoints,4,dimnames=list(NULL,c('estimate',ui_interval,'covariate'))))
      names(fitted_switch) <- paste0('T',1:K)
    }

    est <- lapply(seq_len(K), function(x) matrix(0.0,1,3,dimnames=list('slope',c('estimate',ui_interval))))
    for (j in seq_len(K)){

      sims <- cmat %*% t(simbetas[[j]])

      if (modelframe_classes[i] == 'numeric'){

        sims_switch <- cmat_switch %*% t(simbetas[[j]])

        est[[j]][,] <- c(mean(simbetas[[j]][,cov_i]),quantile(simbetas[[j]][,cov_i],c(ui_offset,1-ui_offset)))

        fitted[[j]][,] <- c(rowMeans(sims),t(apply(sims,1,function(x) quantile(x,c(ui_offset,1-ui_offset)))),cdata[,cov_i])
        fitted_switch[[j]][,] <- c(rowMeans(sims_switch),t(apply(sims_switch,1,function(x) quantile(x,c(ui_offset,1-ui_offset)))),cdata_switch[,cov_i])

      }else if (modelframe_classes[i] == 'factor'){

        diff <- sims[1,] - sims[2,]
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

estimate_topic_effects <- function(covariate,topics,metadata,nsims=100,ci=.95){

  fit <- topics$fit
  K <- fit$settings$dim$K
  formula <- fit$settings$covariates$formula

  if (!is.null(formula)){
    formula <- as.formula(sprintf('1:%s %s',K,paste0(formula,collapse=' ')))
  }else{
    formula <- as.formula(sprintf('1:%s ~ %s',K,covariate))
  }

  effects <- stm::estimateEffect(formula,fit,metadata,uncertainty='Global')

  cov_type <- length(unique(effects$modelframe[,covariate]))

  if (cov_type > 2){
    covariate_effects <- estimate_topic_effects_continuous(effects,covariate,nsims=nsims,ci=ci)
  }else if (cov_type == 2){
    covariate_effects <- estimate_topic_effects_binary(effects,covariate,nsims=nsims,ci=ci)
  }else{
    stop('Covariate has only 1 level!')
  }

}

estimate_topic_effects_binary <- function(effects,covariate,nsims=100,ci=.95){

  K <- length(effects$topics)

  covariate_levels <- levels(as.factor(effects$modelframe[,covariate]))

  cthis <- stm:::produce_cmatrix(prep=effects,covariate=covariate,
                           method='difference',
                           cov.value1=covariate_levels[1],cov.value2=covariate_levels[2],
                           npoints=100)
  cdata <- cthis$cdata
  cmat <- cthis$cmatrix
  simbetas <- stm:::simBetas(effects$parameters,nsims=nsims)
  offset <- (1-ci)/2

  uvals <- levels(cdata[,covariate])

  est <- lapply(seq_len(K),function(x) matrix(0.0,length(uvals),3,
                                              dimnames=list(uvals,c('estimate',paste0(c(offset,1-offset)*100,'%')))))
  for (i in seq_len(K)){

    sims <- cmat %*% t(simbetas[[i]])
    diff <- sims[1,] - sims[2,]
    est[[i]][,1] <- mean(diff)
    est[[i]][,2:3] = quantile(diff,c(offset,1-offset))

  }

  return(est)

}

estimate_topic_effects_continuous <- function(effects,covariate,nsims=100,ci=.95){

  K <- length(effects$topics)

  cthis <- stm:::produce_cmatrix(prep=effects,covariate=covariate,
                                 method='continuous',npoints=100)
  cdata <- cthis$cdata
  cmat <- cthis$cmatrix
  simbetas <- stm:::simBetas(effects$parameters,nsims=nsims)
  offset <- (1-ci)/2

  uvals <- cdata[,covariate]

  est <- matrix(0.0,K,length(uvals))
  cis <- matrix(0.0,K,2,dimnames=list(NULL,paste0(c(offset,1-offset)*100,'%')))

  est <- lapply(seq_len(K),function(x) matrix(0.0,length(uvals),3,
                                              dimnames=list(NULL,c('estimate',paste0(c(offset,1-offset)*100,'%')))))
  for (i in seq_len(K)){

      sims <- cmat %*% t(simbetas[[i]])
      est[[i]][,1] <- rowMeans(sims)
      est[[i]][,2:3] = apply(sims,1,function(x) quantile(x,c(offset,1-offset)))

  }

  return(est)

}

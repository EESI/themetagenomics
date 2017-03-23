#' Estimate topic effects
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

estimate_topic_effects <- function(topics,metadata,formula,refs,nsims=100,ui_level=.8,npoints=100,verbose=FALSE){

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

  formula <- as.formula(sprintf('1:%s %s',K,paste0(formula,collapse=' ')))

  if (verbose) cat('Estimating regression weights with global uncertainty.\n')
  estimated_effects <- stm::estimateEffect(formula,fit,modelframe,uncertainty='Global')

  estimated_effects$modelframe_full <- topics$modelframe # maybe remove htis
  estimated_effects$modelframe <- modelframe

  topic_effects <- estimate_topic_effects_backend(estimated_effects,nsims=nsims,ui_level=ui_level,npoints=npoints,verbose=verbose)

  out <- list(topic_effects=topic_effects,topics=topics,modelframe=topics$modelframe)
  class(out) <- 'effects'
  attr(out,'type') <- 'topics'

  return(out)

}

#' Backend to extract effects for \code{\link{estimate_topic_effects}}.
#' @keywords internal
estimate_topic_effects_backend <- function(estimated_effects,nsims=100,ui_level=.8,npoints=100,verbose=FALSE){

  K <- length(estimated_effects$topics)

  ui_offset <- (1-ui_level)/2
  ui_interval <- paste0(c(ui_offset,1-ui_offset)*100,'%')
  ui_levels <- matrix(0.0,K,2,dimnames=list(NULL,ui_interval))

  modelframe <- estimated_effects$modelframe
  modelframe_full <- estimated_effects$modelframe_full

  if (verbose) cat('Simulating beta coeffiicents from MVN.\n')
  simbetas <- stm:::simBetas(estimated_effects$parameters,nsims=nsims)
  for (i in seq_along(simbetas)) colnames(simbetas[[i]]) <- names(estimated_effects$parameter[[1]][[1]]$est)

  multiclasses <- create_multiclasses_table(modelframe,modelframe_full)
  mods <- multiclasses[multiclasses[,'class'] == 'factor','full']

  covariate_list <- vector(mode='list',length=ncol(estimated_effects$modelframe_full))
  names(covariate_list) <- colnames(estimated_effects$modelframe_full)
  for (i in seq_len(nrow(multiclasses))){

    cov_i <- multiclasses[i,'full']
    if (verbose) cat(sprintf('Making posterior predictions for %s.\n',cov_i))

    if (multiclasses[i,'class'] == 'numeric' & length(mods) > 0){

      est <- matrix(0.0,K,3,dimnames=list(1:K,c('estimate',ui_interval)))
      for (j in seq_len(K)) est[j,] <- c(mean(simbetas[[j]][,cov_i]),quantile(simbetas[[j]][,cov_i],c(ui_offset,1-ui_offset)))

      # Preallocate lists
      fitted <- lapply(seq_len(K),function(x) matrix(0.0,npoints,4,dimnames=list(NULL,c('estimate',ui_interval,'covariate'))))
      names(fitted) <- paste0('T',1:K)
      fitted <- lapply(seq_along(mods),function(x) fitted)
      names(fitted) <- mods
      fitted_switch <- lapply(seq_len(K),function(x) matrix(0.0,npoints,4,dimnames=list(NULL,c('estimate',ui_interval,'covariate'))))
      names(fitted_switch) <- paste0('T',1:K)
      fitted_switch <- lapply(seq_along(mods),function(x) fitted_switch)
      names(fitted_switch) <- mods

      for (mod in mods){

        if (verbose) cat(sprintf('Making posterior predictions for %s given %s.\n',cov_i,mod))

        ppdx <- make_ppd_x(covariate=cov_i,mod=mod,modelframe_full=modelframe_full,npoints=npoints)

        for (j in seq_len(K)){
          sims <- ppdx[[1]][[1]] %*% t(simbetas[[j]])
          sims_switch <- ppdx[[1]][[2]] %*% t(simbetas[[j]])

          fitted[[mod]][[j]][,] <- c(rowMeans(sims),
                                     t(apply(sims,1,function(x) quantile(x,c(ui_offset,1-ui_offset)))),
                                     ppdx[[1]][[1]][,cov_i])
          fitted_switch[[mod]][[j]][,] <- c(rowMeans(sims_switch),
                                            t(apply(sims_switch,1,function(x) quantile(x,c(ui_offset,1-ui_offset)))),
                                            ppdx[[1]][[1]][,cov_i])

        }

      }

    }else{

      fitted <- NULL
      fitted_switch <- NULL

      ppdx <- make_ppd_x(covariate=cov_i,modelframe_full=modelframe_full,npoints=npoints)

      for (j in seq_len(K)){

        sims <- ppdx[[1]] %*% t(simbetas[[j]])
        diff <- sims[2,] - sims[1,]
        est[j,] <- c(mean(diff),quantile(diff,c(ui_offset,1-ui_offset)))

      }

    }

    rownames(est) <- paste0('T',1:K)
    rank <- dense_rank(est[,1])
    names(rank) <- paste0('T',1:K)
    sig <- which(rowSums(sign(est[,2:3])) != 0)

    covariate_list[[i]] <- list(est=est,rank=rank,sig=sig,fitted=fitted,fitted_switch=fitted_switch)

  }

  return(covariate_list)

}

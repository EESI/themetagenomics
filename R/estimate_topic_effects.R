#' Estimate topic effects
#'
#' Given a covariate of interest, measure its relationship with the samples over
#' topics distribution from the STM.
#'
#' @param object (required) Ouput of \code{\link{find_topics}}.
#' @param metadata Matrix or dataframe containing sample information with row or
#'   column names corresponding to the otu_table.
#' @param formula New formula for covariates of interest found in metadata,
#'   different than the formula used to generate object. Interactions,
#'   transformations, splines, and polynomial expansions are permitted.
#' @param refs Character vector of length equal to the number of factors or
#'   binary covariates in formula, indicating the reference level.
#' @param nsims Number of simulations to perform for estimating
#'   covariate effects. Defaults to 100.
#' @param ui_level Width of uncertainty interval for reporting effects. Defaults to
#'   .95.
#' @param npoints Number of posterior predictive samples to draw. Defaults to 100.
#' @param seed Seed for the random number generator to reproduce previous
#'   results.
#' @param verbose Logical flag to print progress information. Defaults to FALSE.
#' @param ... Additional arguments for methods.
#'
#' @return An object of class effects containing
#' \describe{
#' \item{topic_effects}{List of the effect estimates for the covariates in formula.}
#' \item{topics}{Object of class topics containing the original output of find_topics.}
#' \item{modelframe}{Original modelframe.}
#' }
#'
#' @details The posterior predictive estimates are calculated depending on the type of covariate. First, all
#' factors are expanded using dummy variables, setting the reference classes as intercepts. For each topic,
#' the topic frequency over samples is regressed against the expanded design matrix. Covariate weights and the
#' variance-covariance matrix is then calculated, which are used to sample new weights using a multivariate
#' normal distribution.
#'
#' The estimation of a specific covariate effect is performed by calculated y-hat from the posterior predictive
#' distribution by holding all covariates other than the target covariate fixed. This is accomplished by
#' marginalizing over the sample data. This fixed design matrix is then multiplied by the weights simulated
#' from the multivariate normal distribution. For a target binary covariate x (which includes expanded factors),
#' effect estimates are defined as the difference between y-hat when x=1 and y-hat when x=0 is calculated, with the reference
#' covariate designated as 1 (hence negative differences imply a strong effect for the reference class). For
#' continuous covariates, the effect estimates are defined as the regression weight for that covariate of interest.
#' To explore the posterior predictive distribution, y-hat is again calculated, but over a vector of values spanning the
#' range of the continuous covariate, with other covariates held fixed as before. Additional y-hat are then calculated
#' while iteratively setting each binary covariate to 0, to explore their influence on the continuous covariate.
#' Nonlinear covariates (e.g., splines) are treated similarly with respect to y-hat. Their effect estimates, however, are
#' calculated by calculating the Spearman rank correlation coefficient between y-hat and y.
#'
#' For each covariate, the effect estimate is returned. y-hat vectors are returned as well for continuous and nonlinear
#' covariates. All effect estimates are ranked in terms of weight or correlation coefficient. Values not overlapping 0 given
#' a user designed level of uncertainty or returned as "significant."
#'
#' @seealso \code{\link[stm]{estimateEffect}}
#'
#' @references
#' Gelman, A. and Hill, J. (2006). Data Analysis Using Regression and
#' Multilevel/Hierarchical Models. Cambridge University Press; 1 edition.
#'
#' Roberts, M.E., Stewart, B.M., Tingley, D., Lucas, C., Leder-Luis,
#' J., Gadarian, S.K., Albertson, B., & Rand, D.G. (2014). Structural topic
#' models for open-ended survey responses. Am. J. Pol. Sci. 58, 1064–1082.
#'
#' @examples
#' formula <- ~DIAGNOSIS
#' refs <- 'CD'
#'
#' dat <- prepare_data(otu_table=GEVERS$OTU,rows_are_taxa=FALSE,tax_table=GEVERS$TAX,
#'                     metadata=GEVERS$META,formula=formula,refs=refs,
#'                     cn_normalize=TRUE,drop=TRUE)
#'
#' \dontrun{
#' topics <- find_topics(dat,K=15)
#' topic_effects <- est(topics)
#' }
#'
#' @aliases est_topics est.topics
#'
#' @export

est.topics <- function(object,metadata,formula,refs,nsims=100,ui_level=.8,npoints=100,
                       seed=object$seeds$next_seed,verbose=FALSE,...){

  set.seed(check_seed(seed))
  next_seed <- sample.int(.Machine$integer.max,1)

  fit <- object$fit
  K <- fit$settings$dim$K

  if (missing(formula)){
    metadata <- object$metadata
    if (is.null(metadata)) stop('Please provide metadata')
    formula <- fit$settings$covariates$formula
    if (is.null(formula)) stop('Please provide a formula.')
  }
  if (missing(refs)) refs <- object$refs

  # ensure new data is OK
  if (nrow(object$metadata) > nrow(metadata))
    stop(paste('Samples in new data are not found in original fitted data.',
               'New data must contain all of the samples originally fit, with no NA values.',
               'Interpolate if necessary.'))
  new_covs <- labels(terms(extract_spline_info(formula,metadata,remove_only=TRUE)))
  metadata <- metadata[rownames(object$metadata),new_covs,drop=FALSE]
  new_data <- prepare_data(otu_table=object$otu_table,rows_are_taxa=FALSE,tax_table=object$tax_table,
                           metadata=metadata,formula=formula,refs=refs,cn_normalize=FALSE,drop=FALSE,
                           verbose=FALSE)
  metadata <- new_data$metadata
  refs <- new_data$refs

  splines <- check_for_splines(formula,metadata)
  if (splines){
    splineinfo <- extract_spline_info(formula,metadata)
    modelframe <- as.data.frame(unclass(model.frame(splineinfo$formula,metadata)))
    modelframe_full <- create_modelframe(splineinfo$formula,metadata,refs)
  }else{
    splineinfo <- NULL
    modelframe <- as.data.frame(unclass(model.frame(formula,metadata)))
    modelframe_full <- create_modelframe(formula,metadata,refs)
  }
  rownames(modelframe) <- rownames(metadata)

  expanded <- expand_multiclass(metadata=modelframe,refs=refs,verbose=verbose)
  modelframe <- expanded$metadata
  refs_type <- expanded$refs_type

  formula <- as.formula(sprintf('1:%s %s',K,paste0(formula,collapse=' ')))

  if (verbose) cat('Estimating regression weights with global uncertainty.\n')
  estimated_effects <- suppressWarnings(stm::estimateEffect(formula,fit,modelframe,uncertainty='Global'))

  estimated_effects$modelframe_full <- modelframe_full
  estimated_effects$modelframe <- modelframe
  estimated_effects$splines <- splineinfo$info

  topic_effects <- est_topics_backend(estimated_effects,fit$theta,
                                      nsims=nsims,ui_level=ui_level,npoints=npoints,verbose=verbose)

  out <- list(topic_effects=topic_effects,topics=object,modelframe=modelframe_full,seeds=list(seed=seed,next_seed=next_seed))
  class(out) <- 'effects'
  attr(out,'type') <- 'topics'
  attr(out,'refs') <- refs_type

  return(out)

}

#' Backend to extract effects for \code{\link{est.topics}}.
#' @keywords internal
est_topics_backend <- function(estimated_effects,theta,nsims=100,ui_level=.8,npoints=100,verbose=FALSE){

  K <- length(estimated_effects$topics)

  ui_offset <- (1-ui_level)/2
  ui_interval <- paste0(c(ui_offset,1-ui_offset)*100,'%')

  if (verbose) cat('Simulating beta coeffiicents from MVN.\n')
  sim_weights <- ppd_weights(estimated_effects$parameters,nsims=nsims)
  for (i in seq_along(sim_weights)) colnames(sim_weights[[i]]) <- names(estimated_effects$parameter[[1]][[1]]$est)

  splineinfo <- estimated_effects$splines
  if (!is.null(splineinfo)) spline_idx <- sapply(splineinfo,function(x) x$var) else spline_idx <- NULL

  multiclasses <- create_multiclasses_table(estimated_effects$modelframe,estimated_effects$modelframe_full,spline_idx)
  mods <- multiclasses[multiclasses$baseclass == 'factor','full']

  covariate_list <- vector(mode='list',length=ncol(estimated_effects$modelframe_full))
  names(covariate_list) <- colnames(estimated_effects$modelframe_full)
  for (i in seq_len(nrow(multiclasses))){

    fitted <- NULL
    fitted_switch <- NULL
    est <- matrix(0.0,K,3,dimnames=list(1:K,c('estimate',ui_interval)))

    cov_i <- multiclasses$full[i]
    attr(cov_i,'baseclass') <- multiclasses$baseclass[i]
    attr(cov_i,'multiclass') <- multiclasses$multiclass[i]

    if (verbose) cat(sprintf('Making posterior predictions for %s.\n',cov_i))

    if (multiclasses$multiclass[i] == 'spline'){
      # for splines, measure posterior sd to get an idea of just how nonlinear each topic is
      ppd <- make_ppd_x(estimated_effects,npoints=100)[[1]]
      for (k in seq_len(K)){
        ppd_beta <- ppd %*% t(sim_weights[[k]])
        ppd_sd <- apply(ppd_beta,2,sd)
        est[k,] <- c(mean(ppd_sd),quantile(ppd_sd,c(ui_offset,1-ui_offset)))
      }
    }else{
      # for non-spline continuous, measure posterior slopes
      for (j in seq_len(K)) est[j,] <- c(mean(sim_weights[[j]][,cov_i]),quantile(sim_weights[[j]][,cov_i],c(ui_offset,1-ui_offset)))
    }

    if (multiclasses$baseclass[i] == 'factor'){
      # for factors, measure the covariate weight (same as ppd diff)
      for (j in seq_len(K)) est[j,] <- c(mean(sim_weights[[j]][,cov_i]),quantile(sim_weights[[j]][,cov_i],c(ui_offset,1-ui_offset)))
    }

    if (multiclasses$baseclass[i] == 'numeric'){

      # cov_vals <- seq(min(estimated_effects$data[[cov_i]]),max(estimated_effects$data[[cov_i]]),length.out=length(estimated_effects$data[[cov_i]]))
      cov_vals <- estimated_effects$data[[cov_i]]

      if (length(mods) > 0){

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

          ppd <- make_ppd_x(estimated_effects,covariate=cov_i,mod=mod,npoints=100)

          for (k in seq_len(K)){
            sims <- ppd[[2]] %*% t(sim_weights[[k]])
            sims_switch <- ppd[[1]] %*% t(sim_weights[[k]])
            fitted[[mod]][[k]] <- cbind(rowMeans(sims),
                                        t(apply(sims,1,function(x) quantile(x,c(ui_offset,1-ui_offset)))),
                                        cov_vals)
            fitted_switch[[mod]][[k]] <- cbind(rowMeans(sims_switch),
                                               t(apply(sims_switch,1,function(x) quantile(x,c(ui_offset,1-ui_offset)))),
                                               cov_vals)
          }

        }

      }else{

        # Preallocate lists
        fitted <- lapply(seq_len(K),function(x) matrix(0.0,npoints,4,dimnames=list(NULL,c('estimate',ui_interval,'covariate'))))
        names(fitted) <- paste0('T',1:K)
        fitted <- list(fitted)

        ppd <- make_ppd_x(estimated_effects,covariate=cov_i,npoints=100)[[1]]

        for (k in seq_len(K)){
          sims <- ppd %*% t(sim_weights[[k]])

          # was fitted[[mod]][[k]] but there isn't a mod
          fitted[[1]][[k]] <- cbind(rowMeans(sims),
                                      t(apply(sims,1,function(x) quantile(x,c(ui_offset,1-ui_offset)))),
                                      cov_vals)
        }

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

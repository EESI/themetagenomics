#' @import lme4 rstan Rcpp
#' @importFrom stats4 summary
NULL

#' Estimate function effects via HMC
#'
#' Given within topic functional predictions, estimate the effects at a given
#' gene function category level via HMC. The effects correspond to a topic-gene
#' category interaction term after accounting for topic and gene category
#' effects.
#'
#' @param gene_table A gene_table formatted via \code{\link{format_gene_table}}.
#' @param inits List of values for parameter initialization. If omitted, values
#'   are generated via \code{\link{glmer.nb}}
#' @param iters Number of iterations for HMC. Defaults to 1000.
#' @param chains (optional) Number of chains for HMC. Defaults to 1.
#' @param return_fit (optional) Logical to return the rstan data and fit.
#' @param verbose (optional) Logical to print progress.
#'
#' @return A list containing
#'
#' \item{summary}{Rstan output that includes the coefficient mean estimates, uncertainty
#' intervals, standard deviations, effective sample size, Rhat statistics}
#' \item{fit}{Rstan fit (if return_fit=TRUE)}
#' \item{data}{Rstan input data (if return_fit=TRUE)}
#' \item{flagged}{Parameter estimates with Rhat statistics > 1.1}
#' \item{inits}{If parameter estimates are flagged, a list of inits to initialize Rstan
#' with more iterations.}
#' @export
resume <- function(object,...) UseMethod('resume')

#' @export
resume.effects <- function(effects_object,init_type=c('last','orig'),inits,
                           iters,chains=1,return_summary=TRUE,verbose=FALSE){

  if (attr(effects_object,'type') != 'functions')
    stop('Effects object must contain functional infrormation.')

  if (missing(inits)){
    init_type <- match.arg(init_type)
    inits <- effects_object$model$inits[[init_type]]
  }

  if (length(inits) < chains)
    inits <- lapply(seq_len(chains),function(x){
      j <- sample(length(inits),1)
      inits[[j]]
      })

  mm <- resume(effects_object$model$fit,
               stan_dat=effects_object$model$data,
               inits=inits,
               gene_table=effects_object$gene_table,pars=effects_object$model$pars,
               iters=iters,chains=chains,return_summary=return_summary,verbose=verbose)

  out <- list(model=mm,gene_table=effects_object$gene_table)
  class(out) <- 'effects'
  attr(out,'type') <- 'functions'
  attr(out,'method') <- attr(effects_object,'method')

  return(out)

}

#' @export
resume.stanfit <- function(stan_obj,stan_dat,inits,gene_table,
                           pars,iters,chains=1,return_summary=TRUE,verbose=FALSE){

  if (chains > 1){
    if (verbose) cat('Preparing parallelization.\n')
    options_old <- options()

    on.exit(options(options_old),add=TRUE)

    rstan::rstan_options(auto_write=TRUE)
    options(mc.cores=chains)
  }

  fit <- rstan::stan(fit=stan_obj,data=stan_dat,
                     init=inits,
                     pars=c('theta'),include=FALSE,
                     iter=iters,chains=chains,
                     verbose=verbose)

  out <- list()
  out[['pars']] <- pars
  out[['fit']] <- fit
  out[['data']] <- stan_dat
  out[['inits']] <- list(orig=inits,
                         last=apply(fit,2,relist,
                              skeleton=rstan:::create_skeleton(fit@model_pars,fit@par_dims)))
  out[['sampler']] <- rstan::get_sampler_params(fit)

  if (return_summary){
    if (verbose) cat('Extracting summary (this often takes some time).\n')
    out[['summary']] <- extract_stan_summary(fit,stan_dat,pars)
    rhat_pars <- pars[pars != 'yhat']
    rhat <- summary(fit,pars=rhat_pars)[['summary']][,'Rhat'] > 1.1
    rhat_count <- sum(rhat,na.rm=TRUE)
    if (rhat_count > 0){
      warning(sprintf('%s parameters with Rhat > 1.1. Consider more iterations.',rhat_count))
      out[['flagged']] <- names(which(rhat))
    }
  }

  return(out)

}

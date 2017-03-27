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

restart_stan_model <- function(gene_table,fit,iters,chains=1,return_fit=TRUE,verbose=FALSE){

  stan_dat <- list(N=nrow(gene_table),
                   J=length(unique(gene_table$pw)),
                   K=length(unique(gene_table$topic)),
                   I=length(levels(as.factor(gene_table$pw):as.factor(gene_table$topic))),
                   pw_full=gene_table$pw,
                   pw=as.integer(as.factor(gene_table$pw)),
                   topic_full=gene_table$topic,
                   topic=gene_table$topic,
                   pwxtopic_full=as.factor(gene_table$pw):as.factor(gene_table$topic),
                   pwxtopic=as.integer(as.factor(gene_table$pw):as.factor(gene_table$topic)),
                   y=gene_table$count)

  fit <- rstan::stan(fit=fit,
                     pars=c('theta'),include=FALSE,
                     iter=iters,chains=chains,
                     verbose=verbose)

  if (verbose) cat('Extracting summary (this often takes some time).\n')

  summary_pars <- c('mu','phi','b_pw_sigma','b_topic_sigma','b_pwxtopic_sigma','b_pw','b_topic','b_pwxtopic','yhat')
  extract_summary <- vector(mode='list',length=length(summary_pars))
  names(extract_summary) <- summary_pars
  for (i in seq_along(extract_summary)){

    if (grepl('^b\\_[a-z]+$',summary_pars[i])){

      extract_summary_tmp <- summary(fit,pars=summary_pars[i],probs=c(.005,.025,.05,.1,.25,.5,.75,0.9,.95,.975,.995))[['summary']]
      par_name <- gsub('^b\\_','',summary_pars[i])
      lookup_table <- unique(cbind(as.character(stan_dat[[par_name]]),as.character(stan_dat[[paste0(par_name,'_full')]])))
      rownames(extract_summary_tmp) <- lookup_table[,2][order(as.integer(lookup_table[,1],decreasing=TRUE))]

      if (par_name == 'pwxtopic'){
        par_name_tmp <- do.call('rbind',strsplit(rownames(extract_summary_tmp),'\\:'))
        colnames(par_name_tmp) <- unlist(strsplit(par_name,'x'))
        extract_summary_tmp <- cbind(par_name_tmp,as.data.frame(extract_summary_tmp))
      }else{
        par_name_tmp <- matrix(rownames(extract_summary_tmp),ncol=1)
        colnames(par_name_tmp) <- par_name
        extract_summary_tmp <- cbind(par_name_tmp,as.data.frame(extract_summary_tmp))
      }

      extract_summary[[i]] <- extract_summary_tmp

    }else{

      extract_summary[[i]] <- summary(fit,pars=summary_pars[i],probs=c(.005,.025,.05,.1,.25,.5,.75,0.9,.95,.975,.995))[['summary']]

    }

  }

  out <- list(summary=extract_summary)

  if (return_fit){
    out[['fit']] <- fit
    out[['data']] <- stan_dat
  }

  rhat <- summary(fit,pars=c('mu','phi','b_pw','b_topic','b_pwxtopic'))[['summary']][,'Rhat'] > 1.1
  rhat_count <- sum(rhat)
  if (rhat_count > 0){

    warning(sprintf('%s parameters with Rhat > 1.1. Consider more iterations.',rhat_count))

    out[['flagged']] <- names(which(rhat))
    out[['inits']] <- list(mu=summary(fit,pars='mu')[['summary']][,'mean'],
                           phi=summary(fit,pars='phi')[['summary']][,'mean'],
                           b_pw=summary(fit,pars='b_pw')[['summary']][,'mean'],
                           b_topic=summary(fit,pars='b_topic')[['summary']][,'mean'],
                           b_pwxtopic=summary(fit,pars='b_pwxtopic')[['summary']][,'mean'])

  }

  out[['sampler']] <- rstan::get_sampler_params(fit)

  return(out)

}

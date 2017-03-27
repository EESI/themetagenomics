#' @import lme4
NULL

#' Estimate function effects via ML
#'
#' Given within topic functional predictions, estimate the effects at a given
#' gene function category level via ML The effects correspond to a topic-gene
#' category interaction term after accounting for topic and gene category
#' effects.
#'
#' @param gene_table A gene_table formatted via \code{\link{format_gene_table}}.
#' @param inits List of values for parameter initialization. If omitted, values
#'   are generated via \link[lme4]{glmer.nb}
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

fit_ml_model <- function(gene_table,iters=1000,return_fit=FALSE,verbose=FALSE,...){

  if (verbose) cat('Fitting model via ML.\n')

  mm <- glmer.nb(count ~ (1|pw) + (1|topic) + (1|pw:topic),
                 data=gene_table,
                 verbose=verbose,
                 control=glmerControl(calc.derivs=TRUE,
                                      optCtrl=list(maxfun=iters)),
                 ...)


  if (verbose) cat('Extracting summary.\n')

  summary_pars <- c('mu','phi','b_pw_sigma','b_topic_sigma','b_pwxtopic_sigma','b_pw','b_topic','b_pwxtopic','yhat')
  extract_summary <- vector(mode='list',length=length(summary_pars))
  names(extract_summary) <- summary_pars

  extract_summary[['mu']] <- matrix(fixef(mm),ncol=1,dimnames=list('mu','mean'))
  extract_summary[['phi']] <- matrix(getME(mm,'glmer.nb.theta'),ncol=1,dimnames=list('phi','mean'))
  extract_summary[['b_pw_sigma']] <- matrix(c(VarCorr(mm)$`pw`),ncol=1,dimnames=list('b_pw_sigma','mean'))
  extract_summary[['b_topic_sigma']] <- matrix(c(VarCorr(mm)$`topic`),ncol=1,dimnames=list('b_topic_sigma','mean'))
  extract_summary[['b_pwxtopic_sigma']] <- matrix(c(VarCorr(mm)$`pw:topic`),ncol=1,dimnames=list('b_pwxtopic_sigma','mean'))
  extract_summary[['b_pw']] <- data.frame(pw=rownames(ranef(mm)$pw),
                                          mean=ranef(mm)$pw[,1])
  rownames(extract_summary[['b_pw']]) <- extract_summary[['pw']]$pw
  extract_summary[['b_topic']] <- data.frame(topic=1:nrow(ranef(mm)$topic),
                                             mean=ranef(mm)$topic[,1])
  rownames(extract_summary[['b_topic']]) <- extract_summary[['topic']]$topic
  extract_summary[['b_pwxtopic']] <- data.frame(pw=gsub('^(.*)\\:([0-9]+)$','\\1',rownames(ranef(mm)$`pw:topic`)),
                                                topic=gsub('^(.*)\\:([0-9]+)$','\\1',rownames(ranef(mm)$`pw:topic`)),
                                                mean=ranef(mm)$`pw:topic`[,1])
  extract_summary[['yhat']] <- matrix(predict(mm),ncol=1)
  dimnames(extract_summary[['yhat']]) <- list(sprintf('yhat[%s]',1:nrow(extract_summary[['yhat']])),'mean')

  out <- list(summary=extract_summary)

  if (return_fit) out[['fit']] <- mm


  return(out)

}

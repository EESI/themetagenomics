#' @importFrom lme4 glmer.nb glmer
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
#'   are generated via \code{\link{lme4::glmer.nb}}
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

fit_ml_model <- function(gene_table,return_fit=FALSE,verbose=FALSE,...){

  mm <- glmer.nb(count ~ (1|pw) + (1|topic) + (1|pw:topic),
                 data=gene_table,...)

  return(mm)

}

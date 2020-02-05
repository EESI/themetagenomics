#' Extract summary statistics
#'
#' @param object Object of class effects, fit via hmc.
#' @param ... Additional arguments for methods.
#'
#' @export
extract <- function(object,...) UseMethod('extract')

#' @describeIn extract Extract summary statistics from HMC effects object
#'
#' Extracts the summary information in a form conducive with vis methods,
#' specifically in cases when return_summary was set to FALSE.
#'
#' @return A list containing
#' \describe{
#' \item{summary}{Rstan summary of parameters from model.}
#' \item{flagged}{Vector of parameter names with Rhat > 1.1.}
#' }
#'
#' @examples
#' formula <- ~DIAGNOSIS
#' refs <- 'Not IBD'
#'
#' dat <- prepare_data(otu_table=GEVERS$OTU,rows_are_taxa=FALSE,tax_table=GEVERS$TAX,
#'                     metadata=GEVERS$META,formula=formula,refs=refs,
#'                     cn_normalize=TRUE,drop=TRUE)
#'
#' \dontrun{
#' topics <- find_topics(dat,K=15)
#'
#' functions <- predict(topics,reference_path='/references/ko_13_5_precalculated.tab.gz')
#'
#' function_effects <- est(functions,level=3,
#'                         iters=500,method='hmc',
#'                         prior=c('laplace','t','laplace'),
#'                         return_summary=FALSE)
#'
#' function_effects_summary <- extract(function_effects)
#' }
#'
#' @export
extract.effects <- function(object,...){

  if (attr(object,'type') != 'functions')
    stop('Effects object must contain functional infrormation.')

  if (attr(object,'method') != 'hmc')
    stop('ML effects object returns summary automatically.')

  fit <- object$model$fit
  pars <- object$model$pars

  out <- list()
  out[['summary']] <- extract_stan_summary(fit,
                                           object$model$data,
                                           pars)
  rhat_pars <- pars[pars != 'yhat']
  rhat <- summary(fit,pars=rhat_pars)[['summary']][,'Rhat'] > 1.1
  rhat_count <- sum(rhat,na.rm=TRUE)
  if (rhat_count > 0){
    warning(sprintf('%s parameters with Rhat > 1.1. Consider more iterations.',rhat_count))
    out[['flagged']] <- names(which(rhat))
  }

  return(out)

}

#' Backend to extract summary from rstan object
#' @keywords internal
extract_stan_summary <- function(fit,stan_dat,summary_pars){
  extract_summary <- vector(mode='list',length=length(summary_pars))
  names(extract_summary) <- summary_pars
  for (i in seq_along(extract_summary)){

    if (grepl('^b\\_[a-z]+$',summary_pars[i])){

      extract_summary_tmp <- stats4::summary(fit,pars=summary_pars[i],probs=c(.005,.025,.05,.1,.25,.5,.75,0.9,.95,.975,.995))[['summary']]
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

      extract_summary[[i]] <- stats4::summary(fit,pars=summary_pars[i],probs=c(.005,.025,.05,.1,.25,.5,.75,0.9,.95,.975,.995))[['summary']]

    }

  }

  return(extract_summary)
}

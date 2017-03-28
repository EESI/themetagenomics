#' @export
extract <- function(object,...) UseMethod('extract')

#' @export
extract.effects <- function(effects_object,verbose=FALSE){

  if (attr(effects_object,'type') != 'functions')
    stop('Effects object must contain functional infrormation.')

  if (attr(effects_object,'method') != 'hmc')
    stop('ML effects object returns summary automatically.')

  if (verbose) cat('Extracting summary (this often takes some time).\n')

  fit <- effects_object$model$fit
  pars <- effects_object$model$pars

  out <- list()
  out[['summary']] <- extract_stan_summary(fit,
                                           effects_object$model$data,
                                           pars)
  rhat_pars <- pars[pars != 'yhat']
  rhat <- summary(fit,pars=rhat_pars)[['summary']][,'Rhat'] > 1.1
  rhat_count <- sum(rhat)
  if (rhat_count > 0){
    warning(sprintf('%s parameters with Rhat > 1.1. Consider more iterations.',rhat_count))
    out[['flagged']] <- names(which(rhat))
  }

  return(out)

}

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

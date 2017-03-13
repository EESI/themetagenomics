#' Fit a topic model on an OTU table
#'
#' Given an OTU table, this function converts the OTU counts across samples into
#' a document format and then fits a structural topic model by wrapping the stm
#' command from the STM package, see \code{\link{stm}}.
#'
#' @param K A positive integer indicating the number of topics to be estimated.
#' @param otu_table Phyloseq object or OTU table as a data frame or matrix with
#'   OTU IDs as row or column names.
#' @param rows_are_taxa Logical whether OTUs are rows and samples are columns
#'   or vice versa.
#' @param formula (optional) Formula object with no response variable or a
#'   matrix containing topic prevalence covariates. Must provide metadata.
#' @param metadata Matrix or dataframe containing the metadata that corresponds
#'   with the formula.
#' @param sigma_prior (optional) Scalar between 0 and 1. This sets the strength
#'   of regularization towards a diagonalized covariance matrix. Setting the
#'   value above 0 can be useful if topics are becoming too highly correlated.
#' @param model (optional) Prefit STM model object to restart an existing model.
#' @param iters (optional) Maximum number of EM iterations.
#' @param tol (optional) Convergence tolerance.
#' @param batches (optional) Number of groups for memoized inference.
#' @param seed (optional) Seed for the random number generator to reproduce
#'   previous results.
#' @param verbose (optional) Logical indicating whether information should be
#'   printed.
#' @param verbose_n (optional) Integer determining the intervals at which labels
#'   are printed.
#' @param control (optional) List of additional parameters control portions of
#'   the optimization. See details.
#' @return A list containing
#'
#' \item{fit}{The topic model fit}
#' \item{docs}{The documents}
#' \item{vocab}{The vocabulary}
#' @export

find_topics <- function(K,otu_table,rows_are_taxa,control=list(),...){

  if (rows_are_taxa == TRUE){

    otu_table <- t(otu_table)

  }

  user_control_params <- names(control)
  control <- control[user_control_params %in% c('gamma.enet','gamma.ic.k','gamma.ic.k',
                                                'nits','burnin','alpha','eta',
                                                'rp.s','rp.p','rp.d.group.size','SpectralRP','maxV')]

  if (length(control) < length(user_control_params)){
    warning(sprintf('Dropping the following control arguments: %s',
                    paste(user_control_params[!(user_control_params %in% names(control))],collapse=' ')))
  }

  vocab <- colnames(otu_table)
  docs <- lapply(seq_len(nrow(otu_table)), function(i) format_to_docs(otu_table[i,],vocab))
  names(docs) <- rownames(otu_table)

  fit <- stm_wrapper(K=K,docs=docs,vocab=vocab,control=control,...)

  return(list(fit=fit,docs=docs,vocab=vocab))

}

#  An STM wrapper for \code{\link{find_topics}}.

stm_wrapper <- function(K,docs,vocab,metadata=NULL,formula=NULL,sigma_prior=0,
                        model=NULL,iters=500,tol=1e-05,batches=1,seed=NULL,
                        verbose=FALSE,verbose_n=5,control=control){

  if (!missing(metadata)){
    if (!is.null(names(docs)) & !is.null(rownames(metadata))){
      if (!identical(names(docs),rownames(metadata))){
        warning('Sample names in OTU table and metadata are inconsistent!')
      }
    }
  }

  fit <- stm::stm(K=K,documents=docs,vocab=vocab,
                  data=metadata,
                  prevalence=formula,
                  init.type='Spectral',
                  sigma.prior=sigma_prior,
                  model=model,
                  max.em.its=iters,emtol=tol,ngroups=batches,
                  seed=seed,verbose=verbose,reportevery=verbose_n,control=control)

  return(fit)

}

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

find_topics <- function(K,otu_table,rows_are_taxa,formula,metadata,refs,control=list(),...){

  if (!missing(formula)){

    if (missing(metadata)) stop('Must provide metadata if a formula is given.\n')

    err <- try(model.frame(formula,data=metadata,na.action=na.fail),silent=TRUE)
    if (class(err) == 'try-error') stop('NA values in metadata. Please remove (see function prepare_data).')

    # metadata <- as.data.frame(unclass(err))
    # rownames(metadata) <- rownames(err)

    # new
    rnames <- rownames(metadata)
    metadata <- as.data.frame(unclass(metadata))
    rownames(metadata) <- rnames
    # /new

    classes <- sapply(metadata,class)

    if (sum(classes == 'factor') > 0){

      if (missing(refs)){
        warning('References are recommended for factors. Using the first level(s).')
        refs <- unlist(lapply(metadata[,classes == 'factor'],function(x) levels(x)[1]))
      }

      if (sum(classes == 'factor') != length(refs)) stop('A reference is required for each factor.')

      ref_check <- all(sapply(seq_along(refs), function(i) refs[i] %in% lapply(metadata[,classes == 'factor'],levels)[[i]]))
      if (!ref_check) stop('Reference(s) not found in factor(s).')

      j <- 1
      for (i in seq_along(classes)){

        if (classes[i] == 'factor'){
          metadata[,i] <- relevel(as.factor(metadata[,i]),ref=refs[j])
          j <- j+1
        }

      }

      # modelframe <- create_modelframe(formula,refs,metadata)

      # new
      splines <- check_for_splines(formula,metadata)
      if (splines){
        spline_info <- extract_spline_info(formula,metadata)
        modelframe <- create_modelframe(spline_info$formula,refs,metadata)
      }else{
        spline_info <- NULL
        modelframe <- create_modelframe(formula,refs,metadata)
      }
      # /new

    }

  }else{

    formula <- NULL
    metadata <- NULL
    modelframe <- NULL

  }

  if (rows_are_taxa == TRUE) otu_table <- t(otu_table)

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

  fit <- stm_wrapper(K=K,docs=docs,vocab=vocab,formula,metadata=metadata,control=control,...)

  out <- list(fit=fit,docs=docs,vocab=vocab,otu_table=otu_table,
              modelframe=modelframe,spline_info=spline_info$info)
  class(out) <- 'topics'

  return(out)

}

#' Wrapper for \code{\link{stm}}
#' @keywords internal
stm_wrapper <- function(K,docs,vocab,formula=NULL,metadata=NULL,sigma_prior=0,
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

#' Print summary for topics class
#' @export
print.topics <- function(topics_object,...){
  cat(sprintf('A %s object containing a topic model with %s topics, %s samples and %s discrete taxa.\n',
              class(topics_object),
              topics_object$fit$settings$dim$K,
              topics_object$fit$settings$dim$N,
              topics_object$fit$settings$dim$V))
}

#' Predict taxonomic functions
#' @export
predict.topics <- function(topics_object,type=c('function'),...){
  type <- match.arg(type)

  if (type == 'function'){
    predict_functions(topics_object$fit,...)
  }
}

#' Prevent object renaming in class topics
#' @export
`names<-.topics` <- function(object,value){
  warning('topics-class objects cannot be renamed.')
  return(object)
}

#' Prevent attribute renaming in class topics
#' @export
`attributes<-.topics` <- function(object,value){
  warning('topics-class attributes cannot be renamed.')
  return(object)
}

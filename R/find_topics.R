#' @importFrom lda lda.collapsed.gibbs.sampler
NULL

#' Perform topic estimation on a themetadata object
#'
#' Given a themetadata object, this function converts the OTU counts across
#' samples into a document format and then fits a structural topic model by
#' wrapping the \link[stm]{stm} function from package stm.
#'
#' @param themetadata_object (required) Ouput of \code{\link{prepare_data}}.
#' @param K (required) A positive integer indicating the number of topics to be
#'   estimated.
#' @param sigma_prior Scalar between 0 and 1. This sets the strength of
#'   regularization towards a diagonalized covariance matrix. Setting the value
#'   above 0 can be useful if topics are becoming too highly correlated.
#'   Defaults to 0.
#' @param model Prefit STM model object to restart an existing model.
#' @param iters Maximum number of EM iterations. Defaults to 500.
#' @param tol Convergence tolerance. Defaults to 1e-5.
#' @param batches Number of groups for memorized inference. Defaults to 1.
#' @param init_type Type of initialization procedure. Defaults to Spectral
#' @param seed Seed for the random number generator to reproduce previous
#'   results.
#' @param verbose Logical flag to print progress information. Defaults to FALSE.
#' @param verbose_n Integer determining the intervals at which labels are
#'   printed.
#' @param control List of additional parameters control portions of the
#'   optimization. See details.
#'
#' @return An object of class topics containing
#' \describe{
#' \item{fit}{STM object containing topic model fit}
#' \item{docs}{Abundance table in document form of length equal to the number of
#' samples. Each element contains 2-row array, where row 1 contains the the
#' vocabulary index of a given taxon and row 2 contains its abundance in that
#' document}
#' \item{vocab}{Character vector containing vocabulary of taxa IDs, where their
#' position corresponds to the document indexes}
#' \item{otu_table}{Original otu_table}
#' \item{tax_table}{Original tax_table}
#' \item{metadata}{Original metadata}
#' \item{ref}{Original covariate references}
#' \item{modelframe}{Original modelframe}
#' \item{splineinfo}{Original splineinfo}
#' }
#'
#' @details Topics are estimated via \link[stm]{stm} from the stm package. The focus
#' of the themetagenomics pipeline is leveraging both abundance and predicted
#' functional information of 16S rRNA sequencing; hence, the pipeline calls for the
#' use of only "prevalence" information (to use stm terminology). This wrapper
#' therefore removes any options pertaining to "content." If the user is interested
#' in exploring the content component of the STM, then the stm package itself is
#' the ideal place to start. Given that only the prevalence component can be
#' manipulated using find_topics, the following additional parameters can be passed
#' to control as a list (adapted from stm documentation):
#' \describe{
#' \item{gamma.enet}{Scalara between 0 and 1 that controls the degree of L1 and L2
#' regularization, where 0 and 1 correspond to ridge and lasso regression. Defaults
#' to 1.}
#' \item{gamma.ic.k}{Method to select the regularization parameter where 2 corresponds
#' to AIC and log(n) is equivalent to BIC. Defaults to 2.}
#' \item{gamma.maxits}{Maximum number of iterations for estimating prevalence. Defaults
#' to 1000.}
#' \item{nits}{For LDA initialization, the number of Gibbs sampling iterations.
#' Defaults to 50.}
#' \item{burnin}{For LDA initialization, the number of burnin iterations. Defaults to
#' 25.}
#' \item{alpha}{For LDA initialization, the samples over topics distribution
#' hyperparameter.}
#' \item{eta}{For LDA initialization, the topics over words distribution hyperparameter.}
#' \item{rp.s}{For spectral initialization, scalar between 0 and 1 that controls the
#' degree sparsity of random projections. Defaults to .05}
#' \item{rp.p}{For spectral initialization, the dimensionality of random projections.
#' Defaults to 3000.}
#' \item{rp.d.group.size}{For spectral initialization, the block size. Defaults to 2000.}
#' \item{maxV}{For spectral initialization, the maximum number of words used during
#' initialization.}
#' }
#'
#' @seealso \code{\link[glmnet]{glmnet}} \code{\link[stm]{stm}}
#'
#' @references
#' Roberts, M.E., Stewart, B.M., Tingley, D., Lucas, C., Leder-Luis,
#' J., Gadarian, S.K., Albertson, B., & Rand, D.G. (2014). Structural topic
#' models for open-ended survey responses. Am. J. Pol. Sci. 58, 1064–1082.
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
#' }
#'
#' @export
find_topics <- function(themetadata_object,K,sigma_prior=0,model=NULL,iters=500,
                        tol=1e-5,batches=1,init_type=c('Spectral','LDA','Random'),
                        seed=themetadata_object$seed,verbose=FALSE,verbose_n=5,
                        control=list()){

  UseMethod('find_topics')

}

#' @export
find_topics.themetadata <- function(themetadata_object,K,sigma_prior=0,model=NULL,iters=500,
                                    tol=1e-5,batches=1,init_type=c('Spectral','LDA','Random'),
                                    seed=themetadata_object$seeds$next_seed,verbose=FALSE,verbose_n=5,
                                    control=list()){

  set.seed(check_seed(seed))
  mod_seed <- sample.int(.Machine$integer.max,1)
  next_seed <- sample.int(.Machine$integer.max, 1)

  if (missing(K) | K <= 0)
    stop('The number of topic K must be supplied and be a positive non-negative integer.')

  otu_table <- themetadata_object$otu_table
  tax_table <- themetadata_object$tax_table
  formula <- themetadata_object$formula
  metadata <- themetadata_object$metadata
  modelframe <- themetadata_object$modelframe
  splineinfo <- themetadata_object$splineinfo
  refs <- themetadata_object$refs

  init_type <- match.arg(init_type)

  user_control_params <- names(control)
  control <- control[user_control_params %in% c('gamma.enet','gamma.ic.k','gamma.maxits',
                                                'nits','burnin','alpha','eta',
                                                'rp.s','rp.p','rp.d.group.size','maxV')]

  if (length(control) < length(user_control_params)){
    warning(sprintf('Dropping the following control arguments: %s',
                    paste(user_control_params[!(user_control_params %in% names(control))],collapse=' ')))
  }

  vocab <- colnames(otu_table)
  docs <- lapply(seq_len(nrow(otu_table)), function(i) format_to_docs(otu_table[i,],vocab))
  names(docs) <- rownames(otu_table)

  fit <- stm_wrapper(K=as.integer(K),docs=docs,vocab=vocab,formula,metadata=metadata,
                     sigma_prior=sigma_prior,model=model,iters=iters,tol=tol,
                     batches=batches,seed=mod_seed,
                     verbose=verbose,verbose_n=verbose_n,
                     init_type=init_type,control=control)

  out <- list(fit=fit,docs=docs,vocab=vocab,otu_table=otu_table,tax_table=tax_table,metadata=metadata,refs=refs,
              modelframe=modelframe,spline_info=splineinfo$info,seeds=list(seed=seed,mod_seed=mod_seed,next_seed=next_seed))
  class(out) <- 'topics'

  attr(out,'cnn') <- attr(themetadata_object,'cnn')
  attr(out,'refs') <- attr(themetadata_object,'refs')

  return(out)

}

#' Wrapper for stm
#' @keywords internal
stm_wrapper <- function(K,docs,vocab,formula=NULL,metadata=NULL,sigma_prior=0,
                        model=NULL,iters=500,tol=1e-05,batches=1,seed=sample.int(.Machine$integer.max,1),
                        verbose=FALSE,verbose_n=5,init_type=c('Spectral','LDA','Random'),
                        control=control){

  init_type <- match.arg(init_type)

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
                  init.type=init_type,
                  sigma.prior=sigma_prior,
                  model=model,
                  max.em.its=iters,emtol=tol,ngroups=batches,
                  seed=seed,verbose=verbose,reportevery=verbose_n,control=control)

  return(fit)

}

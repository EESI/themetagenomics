#' Estimate the predicted function-topic effects
#'
#' Given within topic functional predictions, estimate the effects at a given
#' gene function category level. The effects correspond to a topic-gene category
#' interaction term after accounting for topic and gene category effects. The
#' model can be fit via either maximum likelihood or Hamiltonian MC -- see
#' details.
#'
#' @param functions List containing the output of
#'   \code{\link{predict_functions}}.
#' @param topics_subset Vector of intergers corresponding to a subset of topics
#'   to be evaluated. Recommended to be < 25.
#' @param level (optional) The gene category level to evalulate. Defaults to 2.
#' @param method (optional) String indicating either ML or HMC. Defaults to ML.
#' @param verbose (optional) Logical to print progress. Defaults to TRUE.
#'
#' @return A glmerMod object (ML) or a list containing rstan output
#'   (HMC).
#' @export

est.functions <- function(functions_object,topics_subset,level=2,method='HMC',...){

  fxn_table <- functions_object$fxn_table
  fxn_meta <- functions_object$fxn_meta

  if (missing(topics_subset) & nrow(fxn_table) <= 25){
    topics_subset <- 1:nrow(fxn_table)
  }

  fxn_table <- fxn_table[topics_subset,]
  fxn_table <- fxn_table[,colSums(fxn_table) > 0]

  fxn_meta <- lapply(fxn_meta,function(x) x[colnames(fxn_table)])

  functions_object$fxn_table <- fxn_table
  functions_object$fxn_meta <- fxn_meta

  gene_table <- format_gene_table(functions_object,level=level)

  if (method == 'ML'){

    mm <- fit_ml_model(gene_table,...)

  }
  if (method == 'HMC'){

    mm <- fit_stan_model(gene_table,...)

  }

  out <- list(model=mm,gene_table=gene_table)
  class(out) <- 'effects'
  attr(out,'type') <- 'functions'
  attr(out,'method') <- method

  return(out)

}

est.effects <- function(effects_object,iters,chains=1,return_fit=TRUE,verbose=FALSE){

  if (attr(effects_object,'type') != 'functions')
    stop('Effects object must contain functional infrormation.')

  if (attr(effects_object,'method') != 'HMC')
    stop('Model must have been fit via HMC.')

  mm <- restart_stan_model(fit=effects_object$model$fit,gene_table=effects_object$gene_table,iters=iters,chains=chains,return_fit=return_fit,verbose=verbose)

  out <- list(model=mm,gene_table=effects_object$gene_table)
  class(out) <- 'effects'
  attr(out,'type') <- 'functions'
  attr(out,'method') <- attr(effects_object,'method')

}

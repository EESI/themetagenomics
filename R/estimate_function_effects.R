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

estimate_function_effects <- function(functions,topics_subset,level=2,method='HML',...){

  fxn_table <- functions$fxn_table
  fxn_meta <- functions$fxn_meta

  if (missing(topics_subset) & nrow(fxn_table) <= 25){
    topics_subset <- 1:nrow(fxn_table)
  }

  fxn_table <- fxn_table[topics_subset,]
  fxn_table <- fxn_table[,colSums(fxn_table) > 0]

  fxn_meta <- lapply(fxn_meta,function(x) x[colnames(fxn_table)])

  functions$fxn_table <- fxn_table
  functions$fxn_meta <- fxn_meta

  gene_table <- format_gene_table(functions,level=level)

  if (method == 'ML'){

    mm <- fit_ml_model(gene_table,...)

  }
  if (method == 'HMC'){

    mm <- fit_stan_model(gene_table,...)

  }

  return(list(model=mm,gene_table=gene_table))

}

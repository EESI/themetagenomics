#' Predict topic functional content
#'
#' Given an STM object, this function predicts the functional content using a
#' precalculated mapping table which maps the taxonomic abundance for a given
#' OTU to functional abundance content across a set of functional genes.
#'
#' @param fit STM object
#' @param reference Location of the precalculated mapping file. See
#'   \code{link{download_ref}}
#' @param scalar (optional) Value for scaling the topics over OTUs distrubution
#'   to predicted counts. Defaults to 100.
#' @param drop (optional) Logical to drop empty gene columns. Defaults to TRUE.
#'
#' @return A list containing
#'
#' \item{fxn_table}{A matrix of gene counts across topics}
#' \item{fxn_meta}{A list of associated functional metadata}
#' \item{pi_meta}{matrix of method specific metadata}
#' @export

predict_functions <- function(fit,reference_path,scalar=100,drop=TRUE,...){

  beta <- round(scalar*exp(fit$beta$logbeta[[1]]))
  rownames(beta) <- paste0('T',1:nrow(beta))
  colnames(beta) <- fit$vocab

  if (grepl('13_5_precalculated',reference_path)){
    predictions <- picrust(beta,rows_are_taxa=FALSE,reference_path=reference_path,drop=drop)
  }else if (grepl('t4f',reference_path)){
    predictions <- t4f(beta,rows_are_taxa=FALSE,reference_path=reference_path,drop=drop,...)
  }else{
    stop('Please provide a valid reference file.')
  }

  return(predictions)

}

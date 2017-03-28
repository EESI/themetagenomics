#' Predict topic functional content
#'
#' Given an object of class topics, this function predicts the functional
#' content using PICRUSt or tax4fun precalculated mapping tables that maps
#' the taxonomic abundance for a given OTU to functional abundance content
#' across a set of functional genes.
#'
#' @param topics_object (required) Output of \code{\link{find_topics}}.
#' @param reference_path Location of the precalculated mapping file, which
#' will determine the method of prediction used.
#' @param scalar Value for scaling the topics over taxa distrubution
#'   to predicted counts. Defaults to 100.
#' @param drop Logical flag to drop empty gene columns. Defaults to TRUE.
#' @param ... Optional arguments for function t4f
#'
#' @return An object of class functions containing
#' \describe{
#' \item{fxn_table}{A matrix of gene counts across topics.}
#' \item{fxn_meta}{A list of functional metadata corresponding to fxn_table.}
#' \item{method_meta}{A matrix of method specific metadata.}
#' }
#'
#' @references
#' ABhauer, K. P., Wemheuer, B. Daniel, R., and Meinicke, P. (2015).
#' Bioinformatics, 1-3. 31(17).
#'
#' Langille, M. G.I.*, Zaneveld, J.*, Caporaso, J. G., McDonald, D., Knights, D.,
#' a Reyes, J., Clemente, J. C., Burkepile, D. E., Vega Thurber, R. L., Knight, R.,
#' Beiko, R. G., and Huttenhower, C. (2013). Nature Biotechnology, 1-10. 8.
#'
#' @seealso \code{\link{download_ref}} \code{\link{picrust}} \code{\link{t4f}}
#'
#' @examples
#' formula <- ~s(age) + drug + sex
#' refs <- c('control','female')
#'
#' dat <- prepare_data(otu_table=OTU,rows_are_taxa=FALSE,tax_table=TAX,
#'                     metadata=META,formula=formula,refs=refs,
#'                     cn_normalize=TRUE,drop=TRUE)
#' topics <- find_topics(dat,K=15)
#'
#' download_ref(destination='/references',reference='gg_ko')
#' functions <- predict(topics,reference_path='/references/ko_13_5_precalculated.tab.gz')
#'
#' @aliases predict_functions predict.functions predict_topics predict.topics
#'
#' @export

predict.topics <- function(topics_object,reference_path,scalar=100,drop=TRUE,...){

  fit <- topics_object$fit

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

  class(predictions) <- 'functions'

  return(predictions)

}

#' Prevent object renaming in class functions
#' @export
`names<-.functions` <- function(object,value){
  warning('functions-class objects cannot be renamed.')
  return(object)
}

#' Prevent attribute renaming in class functions
#' @export
`attributes<-.functions` <- function(object,value){
  warning('functions-class attributes cannot be renamed.')
  return(object)
}


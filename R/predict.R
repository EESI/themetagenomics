#' Predict topic functional content
#'
#' Given an object of class topics, this function predicts the functional
#' content using PICRUSt or tax4fun precalculated mapping tables that maps
#' the taxonomic abundance for a given OTU to functional abundance content
#' across a set of functional genes.
#'
#' @param object (required) Output of \code{\link{find_topics}}.
#' @param reference A string for either gg_ko, gg_cog, or silva_ko.
#' Defaults to gg_ko.
#' @param reference_path Folder path of the reference file
#' @param scalar Value for scaling the topics over taxa distribution
#'   to predicted counts. Defaults to 100.
#' @param cn_normalize Logical flag for performing 16S rRNA copy number
#'   normalization. Defaults to FALSE.
#' @param sample_normalize Logical flag to normalize functional
#'   predictions by the total functional abundance in a sample.
#'   Defaults to FALSE.
#' @param drop Logical flag to drop empty gene columns. Defaults to TRUE.
#' @param ... Additional arguments for t4f method.
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
#' download_ref(destination='/references',reference='gg_ko')
#' functions <- predict(topics,reference='gg_ko',
#'                      reference_path='/references')
#' }
#'
#' @aliases predict_functions predict.functions predict_topics predict.topics
#'
#' @export

predict.topics <- function(object,
                           reference=c('gg_ko','gg_cog','silva_ko'),
                           reference_path,
                           scalar=100,
                           cn_normalize=FALSE,sample_normalize=FALSE,
                           drop=TRUE,...){

  reference <- match.arg(reference)

  if (cn_normalize){
    if (attr(object,'cnn')){
      warning('Copy numbers already normalized via prepare_data. Switching cn_normalize to FALSE.')
      cn_normalize <- FALSE
    }
  }else{
    if (attr(object,'cnn') == FALSE){
      warning('Copy numbers have yet to be normalized. Conisdering cn_normalize=TRUE.')
    }
  }

  fit <- object$fit

  beta <- round(scalar*exp(fit$beta$logbeta[[1]]))
  rownames(beta) <- paste0('T',1:nrow(beta))
  colnames(beta) <- fit$vocab

  if (grepl('gg',reference))
    predictions <- picrust(beta,rows_are_taxa=FALSE,
                           reference=reference,reference_path=reference_path,
                           cn_normalize=cn_normalize,
                           sample_normalize=sample_normalize,
                           drop=drop)

  if (grepl('silva',reference))
    predictions <- t4f(beta,rows_are_taxa=FALSE,
                       tax_table=object$tax_table,
                       reference_path=reference_path,
                       cn_normalize=cn_normalize,
                       sample_normalize=sample_normalize,
                       scalar=scalar,
                       drop=drop,...)

  predictions[['seeds']] <- object$seeds

  class(predictions) <- 'functions'
  if (grepl('ko\\_',reference_path)){
    attr(predictions,'method') <- 'PICRUSt'
    attr(predictions,'ref') <- 'GreenGenes16_5'
    attr(predictions,'db') <- 'KEGG'
  }else if (grepl('cog\\_',reference_path)){
    attr(predictions,'method') <- 'PICRUSt'
    attr(predictions,'ref') <- 'GreenGenes16_5'
    attr(predictions,'db') <- 'COG'
  }else if (grepl('t4f\\_',reference_path)){
    attr(predictions,'method') <- 'Tax4Fun'
    attr(predictions,'ref') <- 'Silva123'
    attr(predictions,'db') <- 'KEGG'
  }else{
    attr(predictions,'type') <- NULL
    attr(predictions,'ref') <- NULL
    attr(predictions,'db') <- NULL
  }
  if (attr(object,'cnn') == FALSE & cn_normalize == FALSE){
    attr(predictions,'cnn') <- FALSE
  }else{
    attr(predictions,'cnn') <- TRUE
  }

  return(predictions)

}


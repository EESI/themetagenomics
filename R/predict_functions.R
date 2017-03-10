#' Predict topic functional content
#'
#' Given an STM object, this function predicts the functional content using
#' a precalculated mapping table which maps the taxonomic abundance for a
#' given OTU to functional abundance content across a set of functional genes.
#'
#' @param fit STM object
#' @param reference Location of the precalculated mapping table. Can be .gz.
#' @param scaler Value for scaling the topics over OTUs distrubution to
#' predicted counts. Defaults at 1000.
#' @return A list containing a matrix of gene counts across topics, a list
#' of associated functional metadata, and matrix of PICRUSt gene metadata.
#' @export

predict_functions <- function(fit,reference_path,scaler=1000){

  beta <- round(scaler*exp(fit$beta$logbeta[[1]]))
  rownames(beta) <- paste0('T',1:nrow(beta))
  colnames(beta) <- fit$vocab

  out <- picrust(reference_path,colnames(beta))
  fxn_mapping <- out$genome_table_out
  rownames(fxn_mapping) <- out$matches

  overlap <- intersect(rownames(fxn_mapping),colnames(beta))
  beta <- beta[,overlap]
  fxn_mapping <- fxn_mapping[overlap,]

  fxn_table <- round(beta %*% fxn_mapping)
  fxn_meta <- format_gene_metadata(out$genemeta)
  pi_meta <- out$pimeta_table_out

  return(list(fxn_table=fxn_table,fxn_meta=fxn_meta,pi_meta=pi_meta))

}

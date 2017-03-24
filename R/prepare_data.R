#' Prepare data for pipeline
#'
#' TBD
#'
#' @param otu_table OTU table.
#' @param rows_are_taxa Logical flag.
#' @param taxa Taxa table.
#' @param metadata Metadata.
#' @param formula Formula for covariates.
#' @param drop Drop empty.
#'
#' @export

prepare_data <- function(otu_table,rows_are_taxa,taxa,metadata,formula,drop=TRUE){

  if (!missing(formula)){
    if (missing(metadata)){
      stop('Must provide metadata if a formula is given.\n')
    }else{
      splines <- check_for_splines(formula,metadata)
      if (splines) formula <- extract_spline_info(formula,metadata,remove_only=TRUE)
      metadata <- model.frame(formula,data=metadata,na.action=na.omit)
    }
  }

  if (rows_are_taxa == TRUE) otu_table <- t(otu_table)

  otu_table <- otu_table[rownames(metadata),]

  if (drop){
    otu_table <- otu_table[,colSums(otu_table) > 0]
    otu_table <- otu_table[rowSums(otu_table) > 0,]
    metadata <- metadata[rownames(otu_table),]
  }

  taxa <- taxa[colnames(OTU),]

  return(list(otu_table=otu_table,taxa=taxa,metadata=metadata))

}

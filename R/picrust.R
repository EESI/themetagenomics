#' Predict topic functional content via PICRUSt
#'
#' Given an OTU abundance table prepared with the GreenGenes reference database,
#' this function predicts the functional content using either COG or KO
#' precalculated mapping tables that map the taxonomic abundance for a given OTU
#' to functional abundance content across a set of functional genes.
#'
#' @param x An abundance table of samples verses OTUs
#' @param rows_are_taxa TRUE/FALSE whether OTUs are rows and samples are columns
#'   or vice versa.
#' @param reference_path Location of the precalculated mapping file. See
#'   \code{link{download_ref}}
#' @param drop (optional) Logical whether to drop empty samples or OTUs after
#'   normalization. Defaults to TRUE.
#'
#' @return A list containing
#'
#' \item{fxn_table}{A matrix of gene counts across topics}
#' \item{fxn_meta}{A list of associated functional metadata}
#' \item{pi_meta}{matrix of PICRUSt metadata (e.g., NSTI)}
#' @export

picrust <- function(x,rows_are_taxa,reference_path,drop=TRUE){

  if (!file.exists(reference_path)){
    stop('Please provide a valid reference file.')
  }

  if (rows_are_taxa == TRUE){

    otu_table <- t(otu_table)

  }

  out <- picrust_otu(reference_path,colnames(x))
  fxn_mapping <- out$genome_table_out
  rownames(fxn_mapping) <- out$matches
  colnames(fxn_mapping) <- out$gene_ids

  overlap <- intersect(rownames(fxn_mapping),colnames(x))
  x <- x[,overlap]
  fxn_mapping <- fxn_mapping[overlap,]

  fxn_table <- round(x %*% fxn_mapping)
  fxn_meta <- format_gene_metadata(out)
  pi_meta <- out$pimeta_table_out
  rownames(pi_meta) <- fit$vocab
  colnames(pi_meta) <- out$pimeta_ids

  if (drop){
    fxn_table <- fxn_table[,colSums(fxn_table)>0]
    fxn_meta <- lapply(fxn_meta,function(x) x[colnames(fxn_table)])
  }

  predictions <- list(fxn_table=fxn_table,fxn_meta=fxn_meta,pi_meta=pi_meta)

  return(predictions)

}

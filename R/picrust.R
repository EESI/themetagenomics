#' @useDynLib themetagenomics, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Predict OTU functional content via PICRUSt
#'
#' Given an OTU abundance table prepared with the GreenGenes reference database,
#' this function predicts the functional content using either COG or KO
#' precalculated mapping tables that map the taxonomic abundance for a given OTU
#' to functional abundance content across a set of functional genes.
#'
#' @param otu_table (required) Matrix or dataframe containing taxa abundances
#'   (counts, non-negative integers) across samples. Rows and columns must be
#'   uniquely named.
#' @param rows_are_taxa (required) Logical flag indicating whether otu_table
#'   rows correspond to taxa (TRUE) or samples (FALSE).
#' @param reference A string for either gg_ko or gg_cog. Defaults to gg_ko.
#' @param reference_path Folder path of the reference file
#' @param cn_normalize Logical flag for performing 16S rRNA copy number
#'   normalization. Defaults to FALSE.
#' @param sample_normalize Logical flag to normalize functional
#'   predictions by the total functional abundance in a sample. Defaults to FALSE.
#' @param drop Logical flag to drop empty gene columns. Defaults to TRUE.
#'
#' @return A list containing
#' \describe{
#' \item{fxn_table}{A matrix of gene counts across topics.}
#' \item{fxn_meta}{A list of functional metadata corresponding to fxn_table.}
#' \item{method_meta}{A matrix of method specific metadata (NSTI).}
#' }
#'
#' @references
#' Langille, M. G.I.*, Zaneveld, J.*, Caporaso, J. G., McDonald, D., Knights, D.,
#' a Reyes, J., Clemente, J. C., Burkepile, D. E., Vega Thurber, R. L., Knight, R.,
#' Beiko, R. G., and Huttenhower, C. (2013). Nature Biotechnology, 1-10. 8.
#'
#' @seealso \code{\link{download_ref}} \code{\link{picrust}}
#'
#' @examples
#' \dontrun{
#' download_ref(destination='/references',reference='gg_ko')
#'
#' predicted_functions <- picrust(otu_table=GEVERS$OTU,rows_are_taxa=TRUE,
#'                                reference='gg_ko',reference_path='/references',
#'                                cn_normalize=TRUE,
#'                                sample_normalize=FALSE,drop=TRUE)
#'                            }
#'
#' @export

picrust <- function(otu_table,rows_are_taxa,reference=c('gg_ko','gg_cog'),reference_path,
                    cn_normalize=FALSE,sample_normalize=FALSE,drop=TRUE){

  reference <- match.arg(reference)

  if (rows_are_taxa == TRUE) otu_table <- t(otu_table)

  if (cn_normalize) otu_table <- cnn(otu_table,rows_are_taxa=FALSE,drop=drop)

  if (reference == 'gg_ko') ref_fn <- 'ko_13_5_precalculated.tab.gz'
  if (reference == 'gg_cog') ref_fn <- 'cog_13_5_precalculated.tab.gz'

  out <- picrust_otu(file.path(reference_path,ref_fn,fsep=platform_sep()),colnames(otu_table))
  fxn_mapping <- out$genome_table_out
  rownames(fxn_mapping) <- out$matches
  colnames(fxn_mapping) <- out$gene_ids

  overlap <- intersect(rownames(fxn_mapping),colnames(otu_table))
  otu_table <- otu_table[,overlap]
  fxn_mapping <- fxn_mapping[overlap,]

  fxn_table <- round(otu_table %*% fxn_mapping)
  fxn_meta <- format_gene_metadata(out)
  pi_meta <- out$pimeta_table_out
  rownames(pi_meta) <- colnames(otu_table)
  colnames(pi_meta) <- out$pimeta_ids

  if (drop){
    fxn_table <- fxn_table[,colSums(fxn_table)>0]
    fxn_meta <- lapply(fxn_meta,function(x) x[colnames(fxn_table)])
  }

  if (sample_normalize) fxn_table <- fxn_table/rowSums(fxn_table)

  predictions <- list(fxn_table=fxn_table,fxn_meta=fxn_meta,method_meta=pi_meta)

  return(predictions)

}

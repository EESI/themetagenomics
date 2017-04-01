#' Normalize an OTU table by 16S rRNA copy number
#'
#' Implements 16S rRNA copy number normalization using the PICRUSt 16S
#' GreenGreenes 13.5 copy number count table (default) or a user provided set of
#' copy numbers.
#'
#' @param otu_table (required) Matrix or dataframe containing taxa abundances
#'   (counts, non-negative integers) across samples. Rows and columns must be
#'   uniquely named.
#' @param rows_are_taxa (required) Logical flag indicating whether otu_table
#'   rows correspond to taxa (TRUE) or samples (FALSE).
#' @param copy_numbers A 2-column matrix or data frame of copy numbers where
#'   column 1 contains the OTU IDs and column 2 the copy numbers.
#' @param drop Logical flag to drop empty rows and columns. Defaults to TRUE.
#' @param verbose Logical flag to print progress information. Defaults to FALSE.
#'
#' @return A normalized, rounded (to nearest integer) abundance table.
#'
#' @references
#' Langille, M. G.I.*, Zaneveld, J.*, Caporaso, J. G., McDonald, D., Knights, D.,
#' a Reyes, J., Clemente, J. C., Burkepile, D. E., Vega Thurber, R. L., Knight, R.,
#' Beiko, R. G., and Huttenhower, C. (2013). Nature Biotechnology, 1-10. 8.
#'
#' @examples
#' nOTU <- cnn(GEVERS$OTU,rows_are_taxa=FALSE,drop=TRUE)
#'
#' @export
cnn <- function(otu_table,rows_are_taxa,copy_numbers,drop=TRUE,verbose=FALSE){

  if (missing(copy_numbers)){
    copy_numbers <- as.matrix(utils::read.table(system.file('references/16S_13_5_precalculated.tab.gz',
                                                     package='themetagenomics'),
                                         sep='\t',header=FALSE,stringsAsFactors=FALSE))
  }else{
    if (ncol(copy_numbers) != 2) stop('copy_numbers must have 2 columns (otu id, copy number).')
    if (diff(colMeans(copy_numbers)) > 0){
      warning('OTU IDs must be in column 1.')
      copy_numbers <- cbind(copy_numbers[,2],copy_numbers[,1])
    }
  }

  if (rows_are_taxa == TRUE) otu_table <- t(otu_table)

  if (is.null(colnames(otu_table))){
    warning('OTU table must have OTU IDs as column names. Returning unnormalized OTU table')
    return(otu_table)
  }
  if (sum(colnames(otu_table) %in% copy_numbers[,1]) == 0){
    warning('OTU IDs must be integer strings of the form GreenGenes 16.X. Returning unnormalized OTU table.')
    return(otu_table)
  }

  otus_to_load <- colnames(otu_table)

  rownames(copy_numbers) <- copy_numbers[,1]
  copy_numbers <- copy_numbers[otus_to_load,]

  norm_table <- round(t(t(otu_table)/copy_numbers[,2]))

  if (drop){
    norm_table <- norm_table[,colSums(norm_table)>0]
    norm_table <- norm_table[rowSums(norm_table)>0,]

    if (verbose){
      cat(sprintf('Dropped %s empty OTU columns.\n',ncol(otu_table)-ncol(norm_table)))
      cat(sprintf('Dropped %s empty sample rows.\n',nrow(otu_table)-nrow(norm_table)))
    }
  }

  return(norm_table)

}

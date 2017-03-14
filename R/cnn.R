#' Normalize an OTU table by 16S rRNA copy number
#'
#' This implements copy number normalization using the PICRUSt 16S GreenGreenes
#' 13.5 copy number count table (default) or a user provided set of copy
#' numbers.
#'
#' @param otu_table Phyloseq object or OTU table as a data frame or matrix with
#'   OTU IDs as row or column names.
#' @param rows_are_taxa TRUE/FALSE whether OTUs are rows and samples are columns
#'   or vice versa.
#' @param copy_numbers A 2-column matrix or data frame of copy numbers where
#'   column 1 contains the OTU IDs and column 2 the copy numbers.
#' @param drop Logical whether to drop empty samples or OTUs after normalization.
#' Defaults to TRUE.
#'
#' @return A normalized, rounded (to nearest integer) OTU table.
#' @export

cnn <- function(otu_table,rows_are_taxa,copy_numbers,drop=TRUE){

  if (missing(copy_numbers)){

    copy_numbers <- as.matrix(read.table(system.file('references/16S_13_5_precalculated.tab.gz',
                                                     package='themetagenomics'),
                                         sep='\t',header=FALSE,stringsAsFactors=FALSE))

  }else{

    if (ncol(copy_numbers) != 2) stop('copy_numbers must have 2 columns (otu id, copy number).')

    if (diff(colMeans(copy_numbers)) > 0){

      warning('otu ids must be in column 1.')
      copy_numbers <- cbind(copy_numbers[,2],copy_numbers[,1])

    }

  }

  if (rows_are_taxa == TRUE){

    otu_table <- t(otu_table)

  }

  if (is.null(rownames(otu_table))) stop('OTU table must have otu ids as row or column names.')

  otus_to_load <- colnames(otu_table)

  rownames(copy_numbers) <- copy_numbers[,1]
  copy_numbers <- copy_numbers[otus_to_load,]

  norm_table <- round(t(t(otu_table)/copy_numbers[,2]))

  if (drop){
    norm_table <- norm_table[,colSums(norm_table)>0]
    norm_table <- norm_table[rowSums(norm_table)>0,]

    cat(sprintf('Dropped %s empty OTU columns.\n',ncol(otu_table)-ncol(norm_table)))
    cat(sprintf('Dropped %s empty sample rows.\n',nrow(otu_table)-nrow(norm_table)))
  }

  return(norm_table)

}

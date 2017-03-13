#' Predict topic functional content via tax4fun
#'
#' Given an abundance table prepared with the Silva reference database, this
#' function predicts the functional content using a KO precalculated mapping
#' table that maps the taxonomic abundance for given taxa to functional
#' abundance content across a set of functional genes.
#'
#' @param x An abundance table of samples verses OTUs
#' @param rows_are_taxa TRUE/FALSE whether OTUs are rows and samples are columns
#'   or vice versa.
#' @param reference_path Location of the precalculated mapping file. See
#'   \code{link{download_ref}}
#' @param type (optional) String for either uproc or pauda. Defaults to uproc.
#' @param short (optional) Logical whether to use a short or long read
#'   reference. Defaults to TRUE.
#' @param copy_number_normalize (optional) Logical whether to perform 16S rRNA
#'   copy number normalization. Defaults to FALSE.
#' @param drop (optional) Logical whether to drop empty samples or OTUs after
#'   normalization. Defaults to TRUE.
#'
#' @return A list containing
#'
#' \item{fxn_table}{A matrix of gene counts across topics}
#' \item{fxn_meta}{A list of associated functional metadata}
#' \item{pi_meta}{matrix of tax4fun metadata (e.g., FTU)}
#' @export

t4f <- function(x,rows_are_taxa,reference_path,type='uproc',short=TRUE,copy_number_normalize=FALSE,drop=TRUE){

  if (rows_are_taxa == TRUE){

    otu_table <- t(otu_table)

  }

  if (short) size <- 'short' else size <- 'long'

  silva_to_kegg <- system.file('extdata/silva_to_kegg.rds',package='themetagenomics')
  copy_numbers <- system.file('extdata/copy_numbers.rds',package='themetagenomics')
  kegg_lookup <- system.file('extdata/kegg_lookup_table.rds',package='themetagenomics')
  ref_profile <- readRDS(reference_path)[[type]][[size]]

  sample_name <- rownames(x)
  otu_ids <- colnames(x)
  silva_ids <- rownames(silva_to_kegg)

  overlap_otus <- intersect(silva_ids,otu_ids)

  sample_sums <- rowSums(otu_table)
  otu_table <- otu_table[,overlap_otus]
  silva_to_kegg <- silva_to_kegg[overlap_otus,]

  tax_to_map <- t(silva_to_kegg) %*% t(otu_table)
  if (copy_number_normalize) tax_to_map <- tax_to_map/copy_numbers
  fxn_table <- t(ref_profile %*% tax_to_map)

  if (drop){
    fxn_table <- fxn_table[,colSums(fxn_table)>0]
    fxn_meta <- lapply(fxn_meta,function(x) x[colnames(fxn_table)])
  }

  desc_tmp <- unique(kegg_lookup_table[,c(1,5)])
  desc_tmp <- desc_tmp[desc_tmp[,1] %in% colnames(fxn_table),]
  desc <- desc_tmp[,1]
  names(desc) <- desc_tmp[,2]

  pw_tmp <- kegg_lookup_table[kegg_lookup_table[,1] %in% colnames(fxn_table),-5]
  pw_kos <- unique(pw_tmp[,1])
  pw <- vector(mode='list',length=length(pw_kos))
  names(pw) <- pw_kos
  for (i in seq_len(nrow(pw_tmp))){
    ko <- pw_tmp[i,1]
    pw[[ko]] <- c(pw[[ko]],list(pw_tmp[i,2:4]))
  }

  fxn_meta <- list(KEGG_Description=desc,KEGG_Pathways=pw)

  ftu <- 1-rowSums(fxn_table)/sample_sums

  return(list(fxn_table=fxn_table,fxn_meta=fxn_meta,ftu=ftu))

  # why do the references have 1 less row than the kegg id
  # why are there blank rows in the kegg id
  # why normalize in the middle
  # why RA normalize
}

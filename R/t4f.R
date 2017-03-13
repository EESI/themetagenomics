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
#' @param taxa Dataframe or matrix containing the taxonomic information (N x 7).
#' @param reference_path Location of the precalculated mapping file. See
#'   \code{link{download_ref}}
#' @param type (optional) String for either uproc or pauda. Defaults to uproc.
#' @param short (optional) Logical whether to use a short or long read
#'   reference. Defaults to TRUE.
#' @param copy_number_normalize (optional) Logical whether to perform 16S rRNA
#'   copy number normalization. Defaults to FALSE.
#' @param sample_normalize (optional) Logical whether to normalize functional
#'   predictions by the total functional abundance in a sample. Defaults to FALSE.
#' @param drop (optional) Logical whether to drop empty samples or OTUs after
#'   normalization and whether to drop KO annotations lacking pathway information.
#'   Defaults to TRUE.
#'
#' @return A list containing
#'
#' \item{fxn_table}{A matrix of gene counts across topics}
#' \item{fxn_meta}{A list of associated functional metadata}
#' \item{pi_meta}{matrix of tax4fun metadata (e.g., FTU)}
#' @export

t4f <- function(taxa_table,rows_are_taxa,taxa,reference_path,type='uproc',short=TRUE,copy_number_normalize=FALSE,sample_normalize=FALSE,drop=TRUE){

  if (missing(taxa)){

    taxa_ids <- colnames(taxa_table)

  }else{

    taxa_ids <- apply(taxa,1,function(l){
      paste0(paste0(l[!is.na(l)],collapse=';'),';')
    })

  }

  if (rows_are_taxa == TRUE){

    taxa_table <- t(taxa_table)

  }

  colnames(taxa_table) <- taxa_ids

  if (short) size <- 'short' else size <- 'long'

  silva_to_kegg <- readRDS(system.file('extdata/silva_to_kegg.rds',package='themetagenomics'))
  copy_numbers <- readRDS(system.file('extdata/copy_numbers.rds',package='themetagenomics'))
  kegg_lookup <- readRDS(system.file('extdata/kegg_lookup_table.rds',package='themetagenomics'))
  ref_profile <- readRDS(reference_path)[[type]][[size]]

  sample_name <- rownames(taxa_table)
  silva_ids <- rownames(silva_to_kegg)

  overlap_taxa <- intersect(silva_ids,taxa_ids)

  sample_sums <- rowSums(taxa_table)
  taxa_table <- taxa_table[,overlap_taxa]
  silva_to_kegg <- silva_to_kegg[overlap_taxa,]

  tax_to_kegg <- t(silva_to_kegg) %*% t(taxa_table)
  if (copy_number_normalize) tax_to_kegg <- tax_to_kegg/copy_numbers
  if (sample_normalize) tax_to_kegg <- t(as.matrix(t(tax_to_kegg))/colSums(tax_to_kegg))
  fxn_table <- t(ref_profile %*% tax_to_kegg)

  if (drop){
    fxn_table <- fxn_table[,colSums(fxn_table)>0]
  }

  desc_tmp <- unique(kegg_lookup[,c(1,5)])
  desc_tmp <- desc_tmp[desc_tmp[,1] %in% colnames(fxn_table),]
  desc <- desc_tmp[,2]
  names(desc) <- desc_tmp[,1]

  if (drop){
    fxn_table <- fxn_table[,names(desc) ]
  }

  pw_tmp <- kegg_lookup[kegg_lookup[,1] %in% colnames(fxn_table),-5]
  pw_kos <- unique(pw_tmp[,1])
  pw <- vector(mode='list',length=length(pw_kos))
  names(pw) <- pw_kos
  for (i in seq_len(nrow(pw_tmp))){
    ko <- pw_tmp[i,1]
    pw[[ko]] <- c(pw[[ko]],list(pw_tmp[i,2:4]))
  }

  fxn_meta <- list(KEGG_Description=desc,KEGG_Pathways=pw)

  t4f_meta <- NULL
  t4f_meta <- cbind(t4f_meta,ftu=1-rowSums(fxn_table)/sample_sums)

  return(list(fxn_table=fxn_table,fxn_meta=fxn_meta,t4f_meta=t4f_meta))

  # why do the references have 1 less row than the kegg id
  # why are there blank rows in the kegg id
  # why RA normalize
}

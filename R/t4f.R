#' Predict taxonomic functional content via tax4fun
#'
#' Given a taxonomic abundance table prepared with the Silva reference database,
#' predicts the functional content using a KO precalculated mapping
#' table that maps the taxonomic abundance for a given tax_table to functional
#' abundance content across a set of functional genes.
#'
#' @param otu_table (required) Matrix or dataframe containing taxa abundances
#'   (counts, non-negative integers) across samples. Rows and columns must be
#'   uniquely named.
#' @param rows_are_taxa (required) Logical flag indicating whether otu_table
#'   rows correspond to taxa (TRUE) or samples (FALSE).
#' @param tax_table Matrix or dataframe containing Silva taxonimic information
#'   with row or column names corresponding to the otu_table. Silva species
#'   information is required.
#' @param reference_path Location of the precalculated mapping file, which
#'   will determine the method of prediction used.
#' @param type Type of protein domain classifcation methods used to generate
#'   references (uproc or pauda). Defaults to uproc.
#' @param short Logical flag whether to use a short or long read
#'   references. Defaults to TRUE.
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
#' \item{method_meta}{A matrix of method specific metadata (FTU).}
#' }
#'
#' @references
#' ABhauer, K. P., Wemheuer, B. Daniel, R., and Meinicke, P. (2015).
#' Bioinformatics, 1-3. 31(17).
#'
#' @seealso \code{\link{download_ref}} \code{\link{picrust}}
#'
#' @examples
#' download_ref(destination='/references',reference='silva_ko')
#'
#' predicted_functions <- t4f(otu_table=OTU,rows_are_taxa=TRUE,reference='/references/t4f_ref_profiles.rds',
#'                            short=TRUE,cn_normalize=TRUE,sample_normalize=TRUE,drop=TRUE)
#'
#' @export

t4f <- function(otu_table,rows_are_taxa,tax_table,reference_path,type='uproc',
                short=TRUE,cn_normalize=FALSE,sample_normalize=FALSE,
                drop=TRUE){

  if (missing(tax_table)){
    taxa_ids <- colnames(otu_table)
  }else{
    taxa_ids <- apply(tax_table,1,function(l){paste0(paste0(l[!is.na(l)],collapse=';'),';')})
  }

  if (rows_are_taxa == TRUE) otu_table <- t(otu_table)

  colnames(otu_table) <- taxa_ids

  if (short) size <- 'short' else size <- 'long'

  silva_to_kegg <- readRDS(system.file('references/silva_to_kegg.rds',package='themetagenomics'))
  copy_numbers <- readRDS(system.file('references/copy_numbers.rds',package='themetagenomics'))
  kegg_lookup <- readRDS(system.file('references/kegg_lookup_table.rds',package='themetagenomics'))
  ref_profile <- readRDS(reference_path)[[type]][[size]]

  sample_name <- rownames(otu_table)
  silva_ids <- rownames(silva_to_kegg)

  overlap_taxa <- intersect(silva_ids,taxa_ids)

  sample_sums <- rowSums(otu_table)
  otu_table <- otu_table[,overlap_taxa]
  silva_to_kegg <- silva_to_kegg[overlap_taxa,]

  tax_to_kegg <- t(silva_to_kegg) %*% t(otu_table)
  if (cn_normalize) tax_to_kegg <- tax_to_kegg/copy_numbers
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

  return(list(fxn_table=fxn_table,fxn_meta=fxn_meta,method_meta=t4f_meta))

  # why do the references have 1 less row than the kegg id
  # why are there blank rows in the kegg id
  # why RA normalize
}

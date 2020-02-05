#' @importFrom Matrix t colSums rowSums
NULL

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
#' @param reference_path Folder path of the silva-to-kegg mapping file
#' (t4f_silva_to_kegg.rds) and reference profiles (t4f_ref_profiles.rds). Must
#' not be renamed.
#' @param type Type of protein domain classification methods used to generate
#'   references (uproc or pauda). Defaults to uproc.
#' @param short Logical flag whether to use a short or long read
#'   references. Defaults to TRUE.
#' @param cn_normalize Logical flag for performing 16S rRNA copy number
#'   normalization. Defaults to FALSE.
#' @param sample_normalize Logical flag to normalize functional
#'   predictions by the total functional abundance in a sample. Defaults to FALSE.
#' @param scalar Value for scaling the topics over functions distrubution
#'   to predicted counts.
#' @param drop Logical flag to drop empty gene columns after prediction. Defaults to TRUE.
#' @param verbose Logical flag to print progress information. Defaults to FALSE.
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
#' \dontrun{
#' download_ref(destination='/references',reference='silva_ko')
#' predicted_functions <- t4f(otu_table=DAVID$OTU,rows_are_taxa=FALSE,
#'                            tax_table=DAVID$TAX,reference='/references',
#'                            type='uproc',short=TRUE,cn_normalize=TRUE,
#'                            sample_normalize=FALSE,scalar=NULL,drop=TRUE)
#'                            }
#'
#' @export

t4f <- function(otu_table,rows_are_taxa,tax_table,reference_path,type=c('uproc','pauda'),
                short=TRUE,cn_normalize=FALSE,sample_normalize=FALSE,scalar,
                drop=TRUE,verbose=FALSE){

  if (any(is.na(otu_table)))
    stop('Please remove NA values from abundance table.')

  type <- match.arg(type)

  ref_profile_path <- file.path(reference_path,'t4f_ref_profiles.rds',fsep=platform_sep())
  map_path <- file.path(reference_path,'t4f_silva_to_kegg.rds',fsep=platform_sep())

  if (rows_are_taxa == TRUE) otu_table <- t(otu_table)

  if (missing(tax_table)){
    taxa_ids <- colnames(otu_table)
  }else{
    tax_table <- tax_table[colnames(otu_table),]
    if (all(grepl('^[a-z]__',tax_table))) tax_table <- gsub('^[a-z]__','',tax_table)
    taxa_ids <- apply(tax_table,1,function(l){paste0(paste0(l[!is.na(l) & l != ''],collapse=';'),';')})
  }
  colnames(otu_table) <- taxa_ids

  if (short) size <- 'short' else size <- 'long'

  silva_to_kegg <- readRDS(map_path)
  copy_numbers <- readRDS(system.file('references/copy_numbers.rds',package='themetagenomics'))
  kegg_lookup <- readRDS(system.file('references/kegg_lookup_table.rds',package='themetagenomics'))
  ref_profile <- readRDS(ref_profile_path)[[type]][[size]]

  silva_ids <- rownames(silva_to_kegg)

  overlap_taxa <- intersect(silva_ids,taxa_ids)

  sample_sums <- rowSums(otu_table)
  otu_table <- otu_table[,overlap_taxa]
  sample_sums_mapped <- rowSums(otu_table)
  silva_to_kegg <- silva_to_kegg[overlap_taxa,]

  tax_to_kegg <- t(silva_to_kegg) %*% t(otu_table)
  if (cn_normalize) tax_to_kegg <- tax_to_kegg/copy_numbers
  if (sample_normalize){
    empty_samps <- colSums(tax_to_kegg) == 0
    if (sum(empty_samps) > 0){
      if (verbose) cat(sprintf('Empty mappings detected, dropping %s samples.\n',sum(empty_samps)))
      tax_to_kegg <- tax_to_kegg[,!empty_samps]
    }
    tax_to_kegg <- t(as.matrix(t(tax_to_kegg))/colSums(tax_to_kegg))
  }
  fxn_table <- t(ref_profile %*% tax_to_kegg)

  desc_tmp <- unique(kegg_lookup[,c(1,5)])
  desc_tmp <- desc_tmp[desc_tmp[,1] %in% colnames(fxn_table),]
  desc <- desc_tmp[,2]
  names(desc) <- desc_tmp[,1]

  if (drop){
    if (verbose) cat('Dropping post-prediction empty rows and columns.\n')
    fxn_table <- fxn_table[,names(desc)]
    fxn_table <- fxn_table[,colSums(fxn_table) > 0]
  }
  t4f_meta <- cbind(NULL,ftu=1-sample_sums_mapped/sample_sums)

  # scale to counts
  if (!missing(scalar)) {
    fxn_table <- round(fxn_table*scalar/max(fxn_table))
    if (drop) fxn_table <- fxn_table[,colSums(fxn_table) > 0]
  }

  pw_tmp <- kegg_lookup[kegg_lookup[,1] %in% colnames(fxn_table),-5]
  pw_kos <- unique(pw_tmp[,1])
  pw <- vector(mode='list',length=length(pw_kos))
  names(pw) <- pw_kos
  for (i in seq_len(nrow(pw_tmp))){
    ko <- pw_tmp[i,1]
    pw[[ko]] <- c(pw[[ko]],list(tolower(pw_tmp[i,2:4])))
  }

  fxn_meta <- list(KEGG_Description=desc,KEGG_Pathways=pw)

  return(list(fxn_table=as.matrix(fxn_table),fxn_meta=fxn_meta,method_meta=t4f_meta))
}

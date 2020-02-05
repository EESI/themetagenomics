#' Inflammatory bowel disease gut microbiome data
#'
#' Subset of samples from the 16S amplicon Gevers et al pediatric
#' inflammatory bowel disease (IBD) dataset. The data include 200 gut samples,
#' 100 of which are controls, spanning 991 OTUs. Three tables are included: an
#' OTU table generated via QIIME, picked against GreenGreens 13.5 at 97%
#' similarity; a taxonomy reference table, and a sample metadata table that
#' includes diagnosis and PCDAI scores, a continuous measure of disease burden.
#'
#' @docType data
#' @name GEVERS
#' @usage GEVERS
#' @keywords datasets
#'
#' @format A list containing a 200x991 matrix (OTU), a 991x7 matrix (TAX), and a
#' 200x3 dataframe (META).
#'
#' @keywords datasets
#'
#' @references Gevers, D., Kugathasan, S., Denson, L.A., et al. (2014). The
#'   Treatment-Naive Microbiome in New-Onset Crohn’s Disease. Cell Host Microbe
#'   15, 382–392. (\href{https://www.ncbi.nlm.nih.gov/pubmed/24629344}{PubMed})
#'
#' @source BioProject: PRJNA237362.
#'   (\href{https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA237362}{PubMed})
#'
#' @examples
#' hist(log(GEVERS$OTU + 1),100)
#' table(GEVERS$META$DIAGNOSIS)
#' boxplot(subset(GEVERS$META,DIAGNOSIS == 'CD')[['PCDAI']])
NULL

#' Human longitudinal microbiome data
#'
#' This includes 16S amplicon sequencing measurements over time from 2
#' individuals. One donor provided both gut and oral samples, whereas the other
#' donor provided only gut samples. The abundance table was generated via Dada2
#' using the Silva reference database. The data span 350 time points.
#'
#' @docType data
#' @name DAVID
#' @usage DAVID
#' @keywords datasets
#'
#' @format A list containing a 746x1493 matrix (ABUND), a 1493x7 matrix (TAX),
#' and a 746x9 dataframe (META).
#'
#' @references David, L. A., Materna, A. C., Friedman, J., Campos-Baptista, M.
#'   I., Blackburn, M. C., Perrotta, A., Erdman, S. E., and Alm, E. J. (2014).
#'   Genome Biology. 15:R89.
#'
#' @source BioProject: PRJEB6518
#'   (\href{https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB6518}{PubMed})
#'
#' @examples
#' hist(log(DAVID$ABUND + 1),100)
#' table(DAVID$META$Site,DAVID$META$Donor)
NULL

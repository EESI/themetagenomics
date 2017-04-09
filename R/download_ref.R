#' Download functional prediction reference tables
#'
#' A function to download the KO and COG 13.5 GreenGenes reference tables for
#' PICRUSt prediction or the KO reference table for tax4fun prediction. The
#' data are stored at \url{https://gitlab.com/sw1/themetagenomics_data/}.
#'
#' @param destination Location of the folder to save the reference files.
#' @param reference A string for either gg_ko, gg_cog, silva_ko, or
#'   all. Defaults to all.
#' @param overwrite Logical flag to overwrite if file already exists. Default
#' to FALSE.
#' @param verbose Logical flag to print progress information. Defaults to FALSE.
#'
#' @seealso \code{\link{picrust}} \code{\link{t4f}}
#'
#' @examples
#' \dontrun{
#' download_ref(destination='/references',reference='gg_ko')
#' }
#'
#' @export

download_ref <- function(destination,reference=c('all','gg_ko','gg_cog','silva_ko'),
                         overwrite=FALSE,verbose=FALSE){

  dir.create(destination,showWarnings=FALSE,recursive=TRUE)

  reference <- match.arg(reference)
  reference_lookup <- c('gg_ko','gg_cog','silva_ko','silva_ko')

  fns <- c('ko_13_5_precalculated.tab.gz',
           'cog_13_5_precalculated.tab.gz',
           't4f_ref_profiles.rds',
           't4f_silva_to_kegg.rds')

  if (reference != 'all')  fns <- fns[reference_lookup %in% reference]

  for (fn in fns){

    if (overwrite==TRUE | !file.exists(file.path(destination,fn))){

      if (verbose) cat(sprintf('Downloading %s.\n',fn))

      suppressWarnings(
      utils::download.file(sprintf('https://gitlab.com/sw1/themetagenomics_data/raw/master/%s',fn),
                    destfile=file.path(destination,fn),
                    method='auto',
                    mode='wb',
                    quiet=!verbose)
      )

    }

  }

}

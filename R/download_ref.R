#' Download PICRUSt reference tables
#'
#' A function to download the KO and COG 13.5 GreenGenes reference tables for
#' PICRUSt prediction.
#'
#' @param reference (optional) A string for either gg_ko, gg_cog, silva_ko, or
#'   all. Defaults to all.
#' @param destination Location of the folder to save the reference files.
#' @param verbose (optional) Logical to print status of download. Defaults to
#'   TRUE.
#'
#' @export

download_ref <- function(destination,reference='all',overwrite=FALSE,verbose=FALSE){

  dir.create(destination,showWarnings=FALSE,recursive=TRUE)

  ref <- c('gg_ko','gg_cog','silva_ko')

  fns <- c('ko_13_5_precalculated.tab.gz',
           'cog_13_5_precalculated.tab.gz',
           't4f_ref_profiles.rds')

  if (reference != 'all'){
    fns <- fns[ref %in% reference]
  }

  for (fn in fns){

    if (overwrite==TRUE | !file.exists(file.path(destination,fn))){

      cat(sprintf('Downloading %s.\n',fn))

      download.file(sprintf('https://github.com/sw1/themetagenomics_data/raw/master/%s',fn),
                    destfile=file.path(destination,fn),
                    quiet=verbose,
                    method='libcurl',
                    mode='w')

    }

  }

}

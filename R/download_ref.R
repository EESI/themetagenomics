#' Download PICRUSt reference tables
#'
#' A function to download the KO and COG 13.5 GreenGenes reference tables
#' for PICRUSt prediction.
#'
#' @param reference A string: ko, cog, or all.
#' @param destination Location of the folder to save the reference files.
#' @param verbose Logical to print status of download.
#' @export

download_ref <- function(reference='all',destination,verbose=FALSE){

  dir.create(destination,showWarnings=FALSE,recursive=TRUE)

  ref <- c('gg_ko','gg_cog','silva_ko')

  fns <- c('ko_13_5_precalculated.tab.gz',
           'cog_13_5_precalculated.tab.gz',
           't4f_ref_profiles.rds')

  if (reference != 'all'){
    fns <- fns[ref %in% reference]
  }

  for (fn in fns){
    download.file(sprintf('https://github.com/sw1/themetagenomics_data/raw/master/%s',fn),
                  destfile(file.path(destination,fn)),
                  quiet=verbose,
                  mode='wb')
  }

}

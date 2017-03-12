#' Download PICRUSt reference tables
#'
#' A function to download the KO and COG 13.5 GreenGenes reference tables
#' for PICRUSt prediction.
#'
#' @param reference A string: ko, cog, or all.
#' @param destination Location of the folder to save the reference files.
#' @param verbose Logical to print status of download.
#' @export

download_ref <- function(reference,destination,verbose=FALSE){

  dir.create(destination,showWarnings=FALSE,recursive=TRUE)

  if (reference == 'all'){
    ref <- c('ko','cog')
  }else if (reference == 'ko'){
    ref <- 'ko'
  }else if (reference == 'cog'){
    ref <- 'cog'
  }

  for (r in ref){
    download.file(sprintf('https://github.com/sw1/themetagenomics_data/raw/master/%s_13_5_precalculated.tab.gz',r),
                  destfile(file.path(destination,sprintf('%s_13_5_precalculated.tab.gz',r))),
                  quiet=verbose,
                  mode='wb')
  }

}

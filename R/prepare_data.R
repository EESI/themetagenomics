#' @importFrom stats model.frame na.omit relevel
NULL

#' Prepare themetadata object from data for topic modeling pipeline
#'
#' Creates a themetadata class by preprocessing data from an OTU table,
#' taxonomic information, sample metadata, and a formula reflecting the preposed
#' relationship between sample metadata and the topics over samples
#' distribution.
#'
#' @param otu_table (required) Matrix or dataframe containing taxa abundances
#'   (counts, non-negative integers) across samples. Rows and columns must be
#'   uniquely named.
#' @param rows_are_taxa (required) Logical flag indicating whether otu_table
#'   rows correspond to taxa (TRUE) or samples (FALSE).
#' @param tax_table Matrix or dataframe containing taxonomic information with row or
#'   column names corresponding to the otu_table.
#' @param metadata Matrix or dataframe containing sample information with row or
#'   column names corresponding to the otu_table.
#' @param formula Formula for covariates of interest found in metadata.
#'   Interactions, transformations, splines, and polynomial expansions are
#'   permitted.
#' @param refs Character vector of length equal to the number of factors or
#'   binary covariates in formula, indicating the reference level.
#' @param cn_normalize Logical flag for performing 16S rRNA copy number
#'   normalization. Defaults to TRUE.
#' @param drop Logical flag to drop empty rows and columns. Defaults to TRUE.
#' @param seed Seed for random number generation. This seed will be passed to
#'   each function that uses this prepared data unless otherwise overridden.
#'   Defaults to a random integer between 1 and the maximum integer supported by R.
#' @param verbose Logical flag to print progress information. Defaults to FALSE.
#'
#' @return An object of class themetadata containing
#' \describe{
#' \item{otu_table}{Matrix of taxa abundances, correctly overlapping with tax_table
#' and metadata. Will be copy number normalized, lacking empty rows and columns by
#' default.}
#' \item{tax_table}{Matrix, correctly overlapping with otu_table}
#' \item{metadata}{Dataframe, correctly overlapping with otu_table and formula. All
#' character covariates are converted to factors.}
#' \item{formula}{Unaltered, given by the user}
#' \item{splineinfo}{List containing the covariate, nonlinear function name, and
#' basis function expansion of all applicable covariates based on the formula.}
#' \item{modelframe}{Dataframe of metadata of only applicable covariates with factors
#' expanded as dummy variables}
#' }
#'
#' @seealso \code{\link[stm]{s}}
#'
#' @examples
#' formula <- ~DIAGNOSIS
#' refs <- 'Not IBD'
#'
#' dat <- prepare_data(otu_table=GEVERS$OTU,rows_are_taxa=FALSE,tax_table=GEVERS$TAX,
#'                     metadata=GEVERS$META,formula=formula,refs=refs,
#'                     cn_normalize=TRUE,drop=TRUE)
#' @export
prepare_data <- function(otu_table,rows_are_taxa,tax_table,metadata,formula,refs,
                         cn_normalize=TRUE,drop=TRUE,seed=sample.int(.Machine$integer.max,1),
                         verbose=FALSE){

  set.seed(check_seed(seed))
  next_seed <- sample.int(.Machine$integer.max,1)

  if (max(otu_table) <= 1)
    stop('Count table must contain counts (non-negative integers) and hence cannot be normalized.')
  if (is.null(colnames(otu_table)) | is.null(rownames(otu_table)))
    stop('otu_table must contain appropriate row and column names.')
  if (!missing(tax_table)){
    if (is.null(colnames(tax_table)) & is.null(rownames(tax_table)))
      stop('tax_table must contain appropriate row and column names.')
    if (!rows_are_taxa){
      if (substr(colnames(otu_table)[1],1,1) == 'X' & substr(colnames(tax_table)[1],1,1) != 'X'){
        warning('otu_table taxa names may have a leading X. Renaming.')
        colnames(otu_table) <- gsub('^X','',colnames(otu_table))
      }
    }
  }
  if (!missing(metadata))
    if (is.null(colnames(metadata)) & is.null(rownames(metadata)))
      stop('metadata must contain appropriate row and column names.')


  slots <- list(otu_table=NULL,
                tax_table=NULL,
                metadata=NULL,
                formula=NULL,
                refs=NULL,
                splineinfo=NULL,
                modelframe=NULL,
                seeds=list(seed=seed,next_seed=next_seed))
  refs_type <- NULL
  splines <- NULL
  class(slots) <- 'themetadata'

  if (rows_are_taxa) {
    otu_table <- t(otu_table)
  }else{
    if (!is.matrix(otu_table)) otu_table <- as.matrix(otu_table)
  }


  miss <- list()
  if (missing(tax_table)) miss$tax_table <- TRUE else miss$tax_table <- FALSE
  if (missing(metadata)) miss$metadata <- TRUE else miss$metadata <- FALSE
  if (missing(formula)) miss$formula <- TRUE else miss$formula <- FALSE
  if (missing(refs) || is.null(refs)){
    miss$refs <- TRUE
    refs <- NULL
  }else{
    miss$refs <- FALSE
  }

  if (!miss$formula)
    if (miss$metadata)
      stop('Must provide metadata if a formula is given.')

  # reorient tax_table and metadata
  if (!miss$metadata){
    if (verbose) cat('Checking row and column names.\n')
    if (sum(rownames(otu_table) %in% colnames(metadata)) > sum(rownames(otu_table) %in% rownames(metadata))){
      if (verbose) cat('Transposing metadata to reorient sample IDs.\n')
      metadata <- t(metadata)
    }
  }
  if (!miss$tax_table & !miss$metadata)
    if (sum(colnames(otu_table) %in% colnames(tax_table)) > sum(colnames(otu_table) %in% rownames(metadata))){
      if (verbose) cat('Transposing tax_table to reorient taxa IDs.\n')
      tax_table <- t(tax_table)
    }

  # make silva like greengenes for downstream functions
  if (!miss$tax_table){
    if (!all(grepl('^[a-z]__',tax_table))){
      tax_table_dimnames <- dimnames(tax_table)
      tax_table[is.na(tax_table)] <- ''
      gg_prefix <- c('k','p','c','o','f','g','s')
      tax_table <- sapply(seq_along(gg_prefix),function(g) paste(gg_prefix[g],tax_table[,g],sep='__'))
      dimnames(tax_table) <- tax_table_dimnames
    }
  }

  # if unsupervised
  if (miss$formula){
    if(is.null(dimnames(otu_table))){
      rownames(otu_table) <- paste0('s',seq_len(nrow(otu_table)))
      colnames(otu_table) <- paste0('f',seq_len(ncol(otu_table)))
      cn_normalize <- FALSE
    }
    if (cn_normalize){
      if (verbose)  cat('Performing copy number normalization.\n')
      otu_table <- cnn(otu_table,FALSE,drop=FALSE)
    }
    if (drop){
      otu_table <- otu_table[,colSums(otu_table) > 0]
      otu_table <- otu_table[rowSums(otu_table) > 0,]
    }
    if (!miss$metadata){
      intersection <- intersect(rownames(otu_table),rownames(metadata))
      otu_table <- otu_table[intersection,]
      metadata <- metadata[intersection,,drop=FALSE]
      slots$metadata <- metadata
    }
    if (!miss$tax_table){
      intersection <- intersect(colnames(otu_table),rownames(tax_table))
      otu_table <- otu_table[,intersection]
      tax_table <- tax_table[intersection,]
      slots$tax_table <- tax_table
    }

    if (verbose) cat('Forcing otu_table to integer mode.\n')
    storage.mode(otu_table) <- 'integer'

    slots$otu_table <- otu_table

    attr(slots,'splines') <- splines
    attr(slots,'refs') <- refs_type
    attr(slots,'cnn') <- cn_normalize
    attr(slots,'drop') <- drop

    return(slots)
  }

  classes <- sapply(metadata,class)
  classes_counts <- c('n'=sum(classes=='numeric' | classes=='integer'),
                      'c'=sum(classes=='character'),
                      'f'=sum(classes=='factor'))
  if (verbose) cat(sprintf('\nStarting stats:
              N otu_table samples: %s
              N otu_table taxa: %s\n
              N metadata numeric %s
              N metadata character %s
              N metadata factor %s\n',
                           nrow(otu_table),
                           ncol(otu_table),
                           classes_counts['n'],
                           classes_counts['c'],
                           classes_counts['f']))
  if (!miss$tax_table)
    if (verbose) cat(sprintf('
              N phyla: %s
              N classes: %s
              N orders: %s
              N families: %s
              N genera: %s
              N species: %s\n\n',
                             length(na.omit(unique(tax_table[,2]))),
                             length(na.omit(unique(tax_table[,3]))),
                             length(na.omit(unique(tax_table[,4]))),
                             length(na.omit(unique(tax_table[,5]))),
                             length(na.omit(unique(tax_table[,6]))),
                             length(na.omit(unique(tax_table[,7])))))

  splines <- check_for_splines(formula,metadata)
  if (splines) formula_tmp <- extract_spline_info(formula,metadata,remove_only=TRUE) else formula_tmp <- formula

  if (verbose) if (any(is.na(metadata))) cat('Removing NA values in metadata.\n')
  metadata <- model.frame(formula_tmp,data=metadata,na.action=na.omit)

  intersection <- intersect(rownames(otu_table),rownames(metadata))
  otu_table <- otu_table[intersection,]
  metadata <- metadata[intersection,,drop=FALSE]

  if (!miss$tax_table){
    intersection <- intersect(colnames(otu_table),rownames(tax_table))
    otu_table <- otu_table[,intersection]
    tax_table <- tax_table[intersection,]
  }
  if (any(is.na(otu_table))) stop('NA values in otu_table. Please correct.\n')

  if (cn_normalize){
    if (verbose)  cat('Performing copy number normalization.\n')
    otu_table <- cnn(otu_table,FALSE,drop=FALSE)
  }

  if (drop){
    otu_table <- otu_table[,colSums(otu_table) > 0]
    otu_table <- otu_table[rowSums(otu_table) > 0,]
    metadata <- metadata[rownames(otu_table),,drop=FALSE]
    if (!miss$tax_table) tax_table <- tax_table[colnames(otu_table),]
  }

  classes <- sapply(metadata,class)
  if(any(classes == 'character')){
    if (verbose) cat('Converting character covariates to factors.\n')
    rnames <- rownames(metadata)
    metadata <- as.data.frame(unclass(metadata))
    rownames(metadata) <- rnames
  }

  expanded <- expand_multiclass(metadata=metadata,refs=refs,verbose=verbose)
  metadata <- expanded$metadata
  refs <- expanded$refs
  refs_type <- expanded$refs_type

  if (splines){
    splineinfo <- extract_spline_info(formula,metadata)
    modelframe <- create_modelframe(splineinfo$formula,metadata,refs)
    slots$splineinfo <- splineinfo
  }else{
    modelframe <- create_modelframe(formula,metadata,refs)
  }
  rownames(metadata) <- rownames(otu_table)

  classes <- sapply(metadata,class)
  classes_counts <- c('n'=sum(classes=='numeric' | classes=='integer'),
                      'c'=sum(classes=='character'),
                      'f'=sum(classes=='factor'))
  if (verbose) cat(sprintf('\nFinal stats:
              N otu_table samples: %s
              N otu_table taxa: %s\n
              N metadata numeric %s
              N metadata character %s
              N metadata factor %s\n',
                           nrow(otu_table),
                           ncol(otu_table),
                           classes_counts['n'],
                           classes_counts['c'],
                           classes_counts['f']))
  if (!miss$tax_table)
    if (verbose) cat(sprintf('
              N phyla: %s
              N classes: %s
              N orders: %s
              N families: %s
              N genera: %s
              N species: %s\n\n',
                             length(na.omit(unique(tax_table[,2]))),
                             length(na.omit(unique(tax_table[,3]))),
                             length(na.omit(unique(tax_table[,4]))),
                             length(na.omit(unique(tax_table[,5]))),
                             length(na.omit(unique(tax_table[,6]))),
                             length(na.omit(unique(tax_table[,7])))))

  if (verbose) cat('Forcing otu_table to integer mode.\n')
  storage.mode(otu_table) <- 'integer'

  slots$otu_table <- otu_table
  if (!miss$tax_table) slots$tax_table <- tax_table
  slots$metadata <- data.frame(metadata,stringsAsFactors=TRUE)
  slots$refs <- refs
  slots$formula <- formula
  slots$modelframe <- modelframe

  attr(slots,'splines') <- splines
  attr(slots,'refs') <- refs_type
  attr(slots,'cnn') <- cn_normalize
  attr(slots,'drop') <- drop

  return(slots)

}

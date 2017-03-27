#' Prepare data for pipeline
#'
#' TBD
#'
#' @param otu_table OTU table.
#' @param rows_are_taxa Logical flag.
#' @param taxa Taxa table.
#' @param metadata Metadata.
#' @param formula Formula for covariates.
#' @param drop Drop empty.
#'
#' @export

prepare_data <- function(otu_table,rows_are_taxa,tax_table,metadata,formula,refs,
                         cn_normalize=TRUE,drop=TRUE,verbose=FALSE){

  slots <- list(otu_table=NULL,
                tax_table=NULL,
                metadata=NULL,
                formula=NULL,
                refs=NULL,
                splineinfo=NULL,
                modelframe=NULL)
  class(slots) <- 'themetadata'

  if (rows_are_taxa) otu_table <- t(otu_table)

  miss <- list()
  if (missing(tax_table)) miss$tax_table <- TRUE else miss$tax_table <- FALSE
  if (missing(metadata)) miss$metadata <- TRUE else miss$metadata <- FALSE
  if (missing(formula)) miss$formula <- TRUE else miss$formula <- FALSE
  if (missing(refs)) miss$refs <- TRUE else miss$refs <- FALSE

  if (!miss$formula)
    if (miss$metadata)
      stop('Must provide metadata if a formula is given.')

  if (miss$formula){
    formula <- ~1
    refs_type <- 'none'
  }
  if (miss$metadata){
    if (rows_raw_taxa){
      metadata <- matrix(1,nrow=ncol(otu_table),2,dimnames=list(colnames(otu_table),NULL))
    }else{
      metadata <- matrix(1,nrow=nrow(otu_table),2,dimnames=list(rownames(otu_table),NULL))
    }
  }
  if (miss$tax_table){
    if (rows_raw_taxa){
      tax_table <- matrix(NA,nrow(otu_table),7,dimnames=list(rownames(otu_table),NULL))
    }else{
      tax_table <- matrix(NA,ncol(otu_table),7,dimnames=list(colnames(otu_table),NULL))
    }
  }

  if (verbose) cat('Checking row and column names.\n')
  if (sum(rownames(otu_table) %in% colnames(metadata)) > sum(rownames(otu_table) %in% rownames(metadata))){
    if (verbose) cat('Transposing metadata to reorient sample IDs.\n')
    metadata <- t(metadata)
  }
  if (sum(colnames(otu_table) %in% colnames(tax_table)) > sum(colnames(otu_table) %in% rownames(metadata))){
    if (verbose) cat('Transposing tax_table to reorient taxa IDs.\n')
    tax_table <- t(tax_table)
  }

  classes <- sapply(metadata,class)
  if (miss$metadata) classes_counts <- c('n'=0,'c'=0,'f'=0) else classes_counts <- c('n'=sum(classes=='numeric'),'c'=sum(classes=='character'),'f'=sum(classes=='factor'))
  if (verbose) cat(sprintf('\nStarting stats:
              N otu_table samples: %s
              N otu_table taxa: %s\n
              N metadata numeric %s
              N metadata character %s
              N metadata factor %s\n
              N phyla: %s
              N classes: %s
              N orders: %s
              N families: %s
              N genera: %s
              N species: %s\n\n',
              nrow(otu_table),
              ncol(otu_table),
              classes_counts['n'],
              classes_counts['c'],
              classes_counts['f'],
              length(na.omit(unique(tax_table[,2]))),
              length(na.omit(unique(tax_table[,3]))),
              length(na.omit(unique(tax_table[,4]))),
              length(na.omit(unique(tax_table[,5]))),
              length(na.omit(unique(tax_table[,6]))),
              length(na.omit(unique(tax_table[,7])))))


  splines <- check_for_splines(formula,metadata)
  if (splines) formula_tmp <- extract_spline_info(formula,metadata,remove_only=TRUE)
  if (!miss$metadata){
    if (verbose) if (any(is.na(metadata))) cat('Removing NA values in metadata.\n')
    metadata <- model.frame(formula_tmp,data=metadata,na.action=na.omit)
  }

  intersection <- intersect(rownames(otu_table),rownames(metadata))
  otu_table <- otu_table[intersection,]
  metadata <- metadata[intersection,]

  intersection <- intersect(colnames(otu_table),rownames(tax_table))
  otu_table <- otu_table[,intersection]
  tax_table <- tax_table[intersection,]

  if (any(is.na(otu_table))) stop('NA values in otu_table. Please correct.\n')
  if (!miss$tax_table) if (any(is.na(tax_table))) stop('NA values in tax_table. Please correct.\n')

  if (cn_normalize){
    if (verbose)  cat('Performing copy number normalization.\n')
    otu_table <- cnn(otu_table,FALSE,drop=FALSE)
  }

  if (drop){
    otu_table <- otu_table[,colSums(otu_table) > 0]
    otu_table <- otu_table[rowSums(otu_table) > 0,]
    metadata <- metadata[rownames(otu_table),]
    tax_table <- tax_table[colnames(otu_table),]
  }

  classes <- sapply(metadata,class)
  if(any(classes == 'character')){
    if (verbose) cat('Converting character covariates to factors.\n')
    metadata <- as.data.frame(unclass(metadata))
  }
  classes <- sapply(metadata,class)
  if (sum(classes == 'factor') > 0){
    if (miss$refs){
      warning('References are recommended for factors. Using the first level(s).')
      refs_type <- 'level_1'
      refs <- unlist(lapply(metadata[,classes == 'factor'],function(x) levels(x)[1]))
    }else{
      refs_type <- 'user'
      if (sum(classes == 'factor') != length(refs)) stop('A reference is required for each factor.')
      ref_check <- all(sapply(seq_along(refs), function(i) refs[i] %in% lapply(metadata[,classes == 'factor'],levels)[[i]]))
      if (!ref_check) stop('Reference(s) not found in factor(s).')

      if (verbose) cat('Setting reference levels for factors.\n')
      j <- 1
      for (i in seq_along(classes)){
        if (classes[i] == 'factor'){
          metadata[,i] <- relevel(as.factor(metadata[,i]),ref=refs[j])
          j <- j+1
        }
      }

      if (splines){
        splineinfo <- extract_spline_info(formula,metadata)
        modelframe <- create_modelframe(splineinfo$formula,refs,metadata)
      }else{
        modelframe <- create_modelframe(formula,refs,metadata)
      }
    }
  }
  rownames(metadata) <- rownames(otu_table)

  classes <- sapply(metadata,class)
  if (miss$metadata) classes_counts <- c('n'=0,'c'=0,'f'=0) else classes_counts <- c('n'=sum(classes=='numeric'),'c'=sum(classes=='character'),'f'=sum(classes=='factor'))
  if (verbose) cat(sprintf('\nFinal stats:
              N otu_table samples: %s
              N otu_table taxa: %s\n
              N metadata numeric %s
              N metadata character %s
              N metadata factor %s\n
              N phyla: %s
              N classes: %s
              N orders: %s
              N families: %s
              N genera: %s
              N species: %s\n\n',
              nrow(otu_table),
              ncol(otu_table),
              classes_counts['n'],
              classes_counts['c'],
              classes_counts['f'],
              length(na.omit(unique(tax_table[,2]))),
              length(na.omit(unique(tax_table[,3]))),
              length(na.omit(unique(tax_table[,4]))),
              length(na.omit(unique(tax_table[,5]))),
              length(na.omit(unique(tax_table[,6]))),
              length(na.omit(unique(tax_table[,7])))))

  slots$otu_table <- otu_table
  slots$modelframe <- modelframe
  slots$splineinfo <- splineinfo
  if (!miss$tax_table) slots$tax_table <- tax_table
  if (!miss$metadata) slots$metadata <- metadata
  if (!miss$refs) slots$refs <- refs
  if (!miss$formula) slots$formula <- formula

  attr(slots,'splines') <- splines
  attr(slots,'refs') <- refs_type
  attr(slots,'cnn') <- cn_normalize
  attr(slots,'drop') <- drop

  return(slots)

}

#' Prevent object renaming in class themetadata
#' @export
`names<-.themetadata` <- function(object,value){
  warning('themetadata-class objects cannot be renamed.')
  return(object)
}

#' Prevent attribute renaming in class themetadata
#' @export
`attributes<-.themetadata` <- function(object,value){
  warning('themetadata-class attributes cannot be renamed.')
  return(object)
}

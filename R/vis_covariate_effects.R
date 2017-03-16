#' Generate interactive binary covariate effect figure
#'
#' TBD
#'
#' @param topics Output of \code{\link{find_topics}} that contains the STM object.
#' @param topic_effects Output of \code{\link{estimate_topic_effects}} that contains regression weights.
#' @param function_effects Output of \code{\link{estimate_function_effects}} that contains the results from either HMC or ML.
#' @param taxa Dataframe or matrix containing the taxonomy information.
#' @param beta_min (optional) Minimum probability in topics over taxa distribution to set to 0.
#' @param gene_min (optional) Mininum count for gene set table.
#'
#' @export


vis_covariate_effects <- function(topics,topic_effects,otu_table,taxa,covariate,metadata,...){

  covariate_data <- metadata[,covariate]

  if (length(unique(covariate_data)) == 2){

    vis_covariate_effects_binary(topics,topic_effects,otu_table,taxa,covariate,metadata,...)

  }else if (length(unique(covariate_data)) > 2){

    if (class(covariate_data) == 'numeric'){

      vis_covariate_effects_continuous(topics,topic_effects,otu_table,taxa,covariate,metadata,...)

    }else{

      # do some multiclass stuff.

    }

  }else{

    stop('Metadata covariate needs to be in the correct form!\n')

  }


}

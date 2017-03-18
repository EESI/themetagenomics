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


vis_covariate_effects <- function(topics,topic_effects,otu_table,taxa,metadata,...){

  if (missing(topic_effects)){

    cat('Estimating topic effects.')

    topic_effects <- estimate_topic_effects(topics,metadata)

  }

  if (sum(is.na(sapply(topic_effects,fitted))) == length(topic_effects)){

    vis_covariate_effects_binary(topics,topic_effects,otu_table,taxa,metadata,...)

  }else if (sum(is.na(sapply(topic_effects,fitted))) == 1){

    vis_covariate_effects_continuous(topics,topic_effects,otu_table,taxa,metadata,...)

  }else{

    stop('Metadata covariate needs to be in the correct form!\n')

  }


}

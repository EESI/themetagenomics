#' Print summary for topics class
#' @export
print.topics <- function(topics_object,...){
  cat(sprintf('A %s object containing a topic model with %s topics, %s samples and %s discrete taxa.\n',
              class(topics_object),
              topics_object$fit$settings$dim$K,
              topics_object$fit$settings$dim$N,
              topics_object$fit$settings$dim$V))
}

#' Print summary for effects class resulting from topics
#' @export
print.effects <- function(effects_object,...){
  cat(sprintf('An %s object containing topic model effects information.\n',class(effects_object)))

  cat(sprintf('\n%s total covariate(s):\n',ncol(effects_object$modelframe)))
  mf <- effects_object$modelframe
  for (i in seq_len(ncol(effects_object$modelframe))){
    cl <- class(mf[[i]])
    if (cl == 'numeric'){
      cat(sprintf('%s\t%s (mean=%.02f, median=%.02f, std=%.02f)\n',
                  cl,
                  colnames(mf)[i],
                  mean(mf[[i]]),median(mf[[i]]),sd(mf[[i]])))
    }else{
      cat(sprintf('%s\t%s (N=%s, N+=%s)\n',
                  cl,
                  colnames(mf)[i],
                  length(mf[[i]]),
                  sum(mf[[i]] == levels(mf[[i]])[2])))
    }
  }

  cat(sprintf('\nTopic weights outside %s-%s uncertainty interval:\n',
              colnames(effects_object$topic_effects[[1]]$est)[2],
              colnames(effects_object$topic_effects[[1]]$est)[3]))
  for (i in seq_len(ncol(effects_object$modelframe))){
    sig_topics <- effects_object$topic_effects[[i]]$sig
    if (length(sig_topics) > 0){
      cat(sprintf('%s:\t%s\n',
                  colnames(mf)[i],
                  paste0('T',sig_topics,collapse=' ')))
    }else{
      cat(sprintf('%s:\t%s\n',
                  colnames(mf)[i],
                  ''))
    }
  }
}

#' Print summary for effects class resulting from functions
#' @export
print.functions <- function(functions_object,...){
  cat(sprintf('Predicted %s functions from %s taxonomic reference database via %s: %s topics across %s genes.\n',
              attr(functions_object,'db'),
              attr(functions_object,'ref'),
              attr(functions_object,'method'),
              nrow(functions_object$fxn_table),
              ncol(functions_object$fxn_table)))
}

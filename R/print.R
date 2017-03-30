#' @importFrom stats median sd
NULL

#' Print summary for topics class
#' @export
print.topics <- function(x,...){
  cat(sprintf('A %s object containing a topic model with %s topics, %s samples and %s discrete taxa.\n',
              class(x),
              x$fit$settings$dim$K,
              x$fit$settings$dim$N,
              x$fit$settings$dim$V))
}

#' Print summary for effects class resulting from topics
#' @export
print.effects <- function(x,...){
  cat(sprintf('An %s object containing topic model effects information.\n',class(x)))

  cat(sprintf('\n%s total covariate(s):\n',ncol(x$modelframe)))
  mf <- x$modelframe
  for (i in seq_len(ncol(x$modelframe))){
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
              colnames(x$topic_effects[[1]]$est)[2],
              colnames(x$topic_effects[[1]]$est)[3]))
  for (i in seq_len(ncol(x$modelframe))){
    sig_topics <- x$topic_effects[[i]]$sig
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
print.functions <- function(x,...){
  cat(sprintf('Predicted %s functions from %s taxonomic reference database via %s: %s topics across %s genes.\n',
              attr(x,'db'),
              attr(x,'ref'),
              attr(x,'method'),
              nrow(x$fxn_table),
              ncol(x$fxn_table)))
}

#' @importFrom stats median sd
NULL

#' @export
print.topics <- function(x,...){
  cat(sprintf('A %s object containing a topic model with %s topics, %s samples and %s discrete taxa.\n',
              class(x)[1],
              x$fit$settings$dim$K,
              x$fit$settings$dim$N,
              x$fit$settings$dim$V))
}

#' @export
print.effects <- function(x,...){

  if (attr(x,'type') == 'topics'){
    cat(sprintf('An %s object containing topic effects information.\n',class(x)[1]))

    cat(sprintf('\n%s total covariate(s):\n',ncol(x$modelframe)))
    mf <- x$modelframe
    for (i in seq_len(ncol(x$modelframe))){
      cl <- class(mf[[i]])[1]
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

  if (attr(x,'type') == 'functions'){

    if (attr(x,'method') == 'hmc'){
      cat(sprintf('An %s object containing function effects information fit via %s.\n',
                  class(x)[1],
                  attr(x,'method')))

      cat(sprintf('%s pathways span %s topics. %s parameters were flagged with Rhat > 1.1 after %s iterations with a %s iteration warmup. The seed for this fit is %s.\n',
                  length(unique(x$gene_table$pw)),
                  length(unique(x$gene_table$topic)),
                  length(x$model$flagged),
                  x$model$fit@stan_args[[1]]$iter,
                  x$model$fit@stan_args[[1]]$warmup,
                  x$model$fit@stan_args[[1]]$seed))

      if (length(x$model$flagged) > 1)
        cat('Given the Rhat results, the model failed to converge. Consider more iterations.\n')
    }

    if (attr(x,'method') == 'ml'){
      cat(sprintf('An %s object containing function effects information fit via %s.\n',
                  class(x)[1],
                  attr(x,'method')))

      cat(sprintf('%s pathways span %s topics. The model was run using %s with %s iterations.\n',
                  length(unique(x$gene_table$pw)),
                  length(unique(x$gene_table$topic)),
                  x$model$fit@optinfo$optimizer,
                  x$model$fit@optinfo$control$maxfun))

      if (grepl('failure',x$model$fit@optinfo$warnings))
        cat('The model failed to converge. Consider more iterations or, preferably, HMC.\n')
    }

  }
}

#' @export
print.functions <- function(x,...){
  cat(sprintf('Predicted %s functions from %s taxonomic reference database via %s: %s topics across %s genes.\n',
              attr(x,'db'),
              attr(x,'ref'),
              attr(x,'method'),
              nrow(x$fxn_table),
              ncol(x$fxn_table)))
}

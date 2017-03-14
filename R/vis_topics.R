#' @importFrom servr httd
NULL

#' Visualize principal components of topics over samples distribution.
#'
#' TBD
#'
#' @param topics Output of \code{\link{find_topics}} that contains the STM object.
#' @param otu_table Taxonomy table used in STM.
#' @param taxa Dataframe or matrix containing the taxonomy information.
#' @param taxa_n (optional) Number of taxa to show in table. Defaults to 25.
#'
#' @export

vis_topics <- function(topics,otu_table,taxa,taxa_n=25){

  K <- topics$fit$settings$dim$K
  vocab <- topics$fit$vocab

  theta <- t(apply(topics$fit$theta + 1/K,1,function(x) x/sum(x)))
  phi <- t(apply(exp(topics$fit$beta$logbeta[[1]]) + .001,1,function(x) x/sum(x)))

  freq <- colSums(otu_table[,vocab])
  len <- sapply(topics$docs,function(x) sum(x[2,]))

  taxon <- pretty_taxa_names(taxa[vocab,])

  json <- createJSON(phi=phi,
                     theta=theta,
                     doc.length=len,
                     vocab=taxon,
                     term.frequency=freq,
                     reorder.topics=FALSE,
                     R=taxa_n)

  tmp <- tempfile()
  suppressMessages(LDAvis::serVis(json,out.dir=tmp,open.browser=TRUE))
  unlink(tmp,recursive=TRUE)

}

# from LDAvis github version
createJSON <- function(phi = matrix(), theta = matrix(), doc.length = integer(),
                       vocab = character(), term.frequency = integer(), R = 30,
                       lambda.step = 0.01, mds.method = LDAvis::jsPCA,
                       cluster,
                       plot.opts = list(xlab = "PC1", ylab = "PC2"),
                       reorder.topics = TRUE,
                       ...) {
  # Set the values of a few summary statistics of the corpus and model:
  dp <- dim(phi)  # should be K x W
  dt <- dim(theta)  # should be D x K

  N <- sum(doc.length)  # number of tokens in the data
  W <- length(vocab)  # number of terms in the vocab
  D <- length(doc.length)  # number of documents in the data
  K <- dt[2]  # number of topics in the model

  # check that certain input dimensions match
  if (dp[1] != K) stop("Number of rows of phi does not match
                       number of columns of theta; both should be equal to the number of topics
                       in the model.")
  if (D != dt[1]) stop("Length of doc.length not equal
                       to the number of rows in theta; both should be equal to the number of
                       documents in the data.")
  if (dp[2] != W) stop("Number of terms in vocabulary does
                       not match the number of columns of phi (where each row of phi is a
                       probability distribution of terms for a given topic).")
  if (length(term.frequency) != W) stop("Length of term.frequency
                                        not equal to the number of terms in the vocabulary.")
  if (any(nchar(vocab) == 0)) stop("One or more terms in the vocabulary
                                   has zero characters -- all terms must have at least one character.")

  # check that conditional distributions are normalized:
  phi.test <- all.equal(rowSums(phi), rep(1, K), check.attributes = FALSE)
  theta.test <- all.equal(rowSums(theta), rep(1, dt[1]),
                          check.attributes = FALSE)
  if (!isTRUE(phi.test)) stop("Rows of phi don't all sum to 1.")
  if (!isTRUE(theta.test)) stop("Rows of theta don't all sum to 1.")

  # compute counts of tokens across K topics (length-K vector):
  # (this determines the areas of the default topic circles when no term is
  # highlighted)
  topic.frequency <- colSums(theta * doc.length)
  topic.proportion <- topic.frequency/sum(topic.frequency)

  # re-order the K topics in order of decreasing proportion:
  if(reorder.topics)
    o <- order(topic.proportion, decreasing = TRUE)
  else
    o <- seq_along(topic.proportion)

  phi <- phi[o, ]
  theta <- theta[, o]
  topic.frequency <- topic.frequency[o]
  topic.proportion <- topic.proportion[o]

  # compute intertopic distances using the specified multidimensional
  # scaling method:
  mds.res <- mds.method(phi)
  if (is.matrix(mds.res)) {
    colnames(mds.res) <- c("x", "y")
  } else if (is.data.frame(mds.res)) {
    names(mds.res) <- c("x", "y")
  } else {
    warning("Result of mds.method should be a matrix or data.frame.")
  }
  mds.df <- data.frame(mds.res, topics = seq_len(K), Freq = topic.proportion*100,
                       cluster = 1, stringsAsFactors = FALSE)
  # note: cluster (should?) be deprecated soon.

  # token counts for each term-topic combination (widths of red bars)
  term.topic.frequency <- phi * topic.frequency

  # compute term frequencies as column sums of term.topic.frequency
  # we actually won't use the user-supplied term.frequency vector.
  # the term frequencies won't match the user-supplied frequencies exactly
  # this is a work-around to solve the bug described in Issue #32 on github:
  # https://github.com/cpsievert/LDAvis/issues/32
  term.frequency <- colSums(term.topic.frequency)
  stopifnot(all(term.frequency > 0))

  # marginal distribution over terms (width of blue bars)
  term.proportion <- term.frequency/sum(term.frequency)

  # Old code to adjust term frequencies. Deprecated for now
  # adjust to match term frequencies exactly (get rid of rounding error)
  #err <- as.numeric(term.frequency/colSums(term.topic.frequency))
  # http://stackoverflow.com/questions/3643555/multiply-rows-of-matrix-by-vector
  #term.topic.frequency <- sweep(term.topic.frequency, MARGIN=2, err, `*`)

  # Most operations on phi after this point are across topics
  # R has better facilities for column-wise operations
  phi <- t(phi)

  # compute the distinctiveness and saliency of the terms:
  # this determines the R terms that are displayed when no topic is selected
  topic.given.term <- phi/rowSums(phi)  # (W x K)
  kernel <- topic.given.term * log(sweep(topic.given.term, MARGIN=2,
                                         topic.proportion, `/`))
  distinctiveness <- rowSums(kernel)
  saliency <- term.proportion * distinctiveness

  # Order the terms for the "default" view by decreasing saliency:
  default.terms <- vocab[order(saliency, decreasing = TRUE)][1:R]
  counts <- as.integer(term.frequency[match(default.terms, vocab)])
  Rs <- rev(seq_len(R))
  default <- data.frame(Term = default.terms, logprob = Rs, loglift = Rs,
                        Freq = counts, Total = counts, Category = "Default",
                        stringsAsFactors = FALSE)
  topic_seq <- rep(seq_len(K), each = R)
  category <- paste0("Topic", topic_seq)
  lift <- phi/term.proportion

  # Collect R most relevant terms for each topic/lambda combination
  # Note that relevance is re-computed in the browser, so we only need
  # to send each possible term/topic combination to the browser
  find_relevance <- function(i) {
    relevance <- i*log(phi) + (1 - i)*log(lift)
    idx <- apply(relevance, 2,
                 function(x) order(x, decreasing = TRUE)[seq_len(R)])
    # for matrices, we pick out elements by their row/column index
    indices <- cbind(c(idx), topic_seq)
    data.frame(Term = vocab[idx], Category = category,
               logprob = round(log(phi[indices]), 4),
               loglift = round(log(lift[indices]), 4),
               stringsAsFactors = FALSE)
  }
  lambda.seq <- seq(0, 1, by=lambda.step)
  if (missing(cluster)) {
    tinfo <- lapply(as.list(lambda.seq), find_relevance)
  } else {
    tinfo <- parallel::parLapply(cluster, as.list(lambda.seq), find_relevance)
  }
  tinfo <- unique(do.call("rbind", tinfo))
  tinfo$Total <- term.frequency[match(tinfo$Term, vocab)]
  rownames(term.topic.frequency) <- paste0("Topic", seq_len(K))
  colnames(term.topic.frequency) <- vocab
  tinfo$Freq <- term.topic.frequency[as.matrix(tinfo[c("Category", "Term")])]
  tinfo <- rbind(default, tinfo)

  # last, to compute the areas of the circles when a term is highlighted
  # we must gather all unique terms that could show up (for every combination
  # of topic and value of lambda) and compute its distribution over topics.

  # unique terms across all topics and all values of lambda
  ut <- sort(unique(tinfo$Term))
  # indices of unique terms in the vocab
  m <- sort(match(ut, vocab))
  # term-topic frequency table
  tmp <- term.topic.frequency[, m]

  # round down infrequent term occurrences so that we can send sparse
  # data to the browser:
  r <- row(tmp)[tmp >= 0.5]
  c <- col(tmp)[tmp >= 0.5]
  dd <- data.frame(Term = vocab[m][c], Topic = r, Freq = round(tmp[cbind(r, c)]),
                   stringsAsFactors = FALSE)

  # Normalize token frequencies:
  dd[, "Freq"] <- dd[, "Freq"]/term.frequency[match(dd[, "Term"], vocab)]
  token.table <- dd[order(dd[, 1], dd[, 2]), ]

  RJSONIO::toJSON(list(mdsDat = mds.df, tinfo = tinfo,
                       token.table = token.table, R = R,
                       lambda.step = lambda.step,
                       plot.opts = plot.opts,
                       topic.order = o))
}

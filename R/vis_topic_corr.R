#' Visualize topic correlations.
#'
#' TBD
#'
#' @param topics Output of \code{\link{find_topics}} that contains the STM object.
#' @param topics_effects Output of \code{\link{estimate_topic_effects}}
#' @param method Correlation graph method. Defaults to huge.
#'
#' @export


vis_topic_corr <- function(topics,topics_effects,method='huge'){

  corr <- stm::topicCorr(topics$fit,method=method)

  effects_rank <- topics_effects$rank
  effects_sig <- topics_effects$sig

  K <- nrow(corr$posadj)
  g <- igraph::graph.adjacency(x$posadj,mode='undirected',
                               weighted=TRUE,diag=FALSE)

  wc <- igraph::cluster_walktrap(g)
  members <- igraph::membership(wc)

  g_d3 <- networkD3::igraph_to_networkD3(g,group=members)

  g_d3$nodes$group <- effects_rank
  g_d3$nodes$sig <- ifelse(1:K %in% effects_sig,50,1)
  g_d3$links <- g_d3$links-1
  g_d3$nodes$name <- paste0('T',g_d3$nodes$name)

  col_scale <- sprintf("color=d3.scaleLinear()\n.domain([1,%s])\n.range(['blue','red']);",K)

  networkD3::forceNetwork(Links=g_d3$links,Nodes=g_d3$nodes,
                          Source='source',Target='target',
                          NodeID='name',Group='group',
                          Nodesize='sig',
                          fontSize=20,
                          zoom=TRUE,
                          colourScale=JS(col_scale))

  return(corr)

}

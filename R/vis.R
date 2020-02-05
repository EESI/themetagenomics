#' @import ggplot2 shiny
#' @importFrom plotly ggplotly plotly plot_ly layout event_data add_trace add_lines add_annotations add_markers plotlyOutput renderPlotly
#' @importFrom stats as.dendrogram dist hclust order.dendrogram cmdscale quantile terms
#' @importFrom scales comma
NULL

#' Launch in interactive visualize to explore topic effects
#'
#' @param object (required) Output of \code{\link{find_topics}}, or
#' \code{\link{est.topics}}, or \code{\link{est.functions}}.
#' @param ... Additional arguments for methods.
#'
#' @details
#'
#' \subsection{Taxa}{
#'   Integrates the samples over topics p(s|k) and topics
#'   over taxa p(k|t) distributions from the STM, the topic correlations from the
#'   p(s|k) component, the covariate effects from the p(s|k) component, and
#'   their relationship with the raw taxonomic abundances. The covariate effects
#'   for each topic are shown as a scatterplot of posterior weights with error bars corresponding the
#'   global approximation of uncertainty. If the covariate chosen is binary, then the posterior regression
#'   weights with uncertainty intervals are shown. This is analogous to the mean difference between factor
#'   levels in the posterior predictive distribution. For continuous covariates, the points
#'   again represent the mean regression weights (i.e., the posterior slope estimate of the
#'   covariate). If, however, a spline or polynomial expansion was used, then the figure shows the posterior estimates of
#'   the standard deviation of the predicted topic probabilities from the posterior predictive distribution.
#'   Colors indicate whether a given point was positive (red) or negative
#'   (blue) and did not enclose 0 at a user defined uncertainty interval.
#'
#'   The ordination figure maintains the color coding just described. The
#'   ordination is performed on p(k|t) via either PCoA (using either
#'   Jensen-Shannon, Euclidean, Hellinger, Bray-Curtis, Jaccard, or Chi-squared
#'   distance) or t-SNE. The latter iterates through decreasing perplexity
#'   values (starting at 30) until the algorithm succeeds. The top 2 or 3
#'   axes can be shown. The radius of the topic points corresponds to the topic
#'   frequencies marginalized over taxa.
#'
#'   The bar plot behaves in accordance with LDAvis. When no topics are chosen,
#'   the overall taxa frequencies are shown. These frequencies do not equal the
#'   abundances found in the initial abundance table. Instead, they show p(k|t)
#'   multiplied by the marginal topic distribution (in counts). To determine the
#'   initial order in which taxa are shown, these two distributions are compared
#'   via Kullback-Liebler divergence and then weighted by the overall taxa
#'   frequency. The coloration of the bars indicates the taxonomic group the
#'   individual taxa belong to. The groups shown are determined based on the
#'   abundance of that group in the raw abundance table. When a topic is
#'   selected, the relative frequency of a given taxa in that topic is shown in
#'   red.
#'
#'   \eqn{\lambda} controls relevance of taxa within a topic, which in turn is used to
#'   adjust the order in which the taxa are shown when a topic is selected.
#'   Relevance is essentially a weighted sum between the probability of taxa in
#'   a given topic and the probability of taxa in a given topic relative to the
#'   overall frequency of that taxa. Adjusting \eqn{\lambda} influences the relative weighting such
#'   that \deqn{r = \lambda x log p(t|k) + \lambda x log p(t|k)/p(x)}
#'   The correlation graph shows the topic correlations from \eqn{p(s|k) ~ MVN(mu,sigma)}.
#'   Again, the coloration described above is conserved. The size
#'   of the nodes reflects the magnitude of the covariate posterior regression weight,
#'   whereas the width of the edges represents the value of the positive
#'   correlation between the connected nodes. By default, the graph estimates
#'   are determined using the the huge package, which first performs a
#'   nonparanormal transformation of p(s|k), followed by a Meinhuasen and
#'   Buhlman procedure. Alternatively, by choosing the simple method, the
#'   correlations are simply a thresholded MAP estimate of p(s|k).
#'   }
#'
#' \subsection{Binary}{
#'   Integrates the topics
#'   over taxa p(k|t) distribution from the STM, binary covariate effects from the p(s|k) component, and
#'   their relationship with the raw taxonomic abundances. The covariate effects
#'   for each topic are shown as a scatterplot of posterior weights with error bars corresponding the
#'   global approximation of uncertainty. If the covariate chosen is binary, then the posterior regression
#'   weights with uncertainty intervals are shown. This is analogous to the mean difference between factor
#'   levels in the posterior predictive distribution. For continuous covariates, the points
#'   again represent the mean regression weights (i.e., the posterior slope estimate of the
#'   covariate). Colors indicate whether a given point was positive (red) or negative
#'   (blue) and did not enclose 0 at a user defined uncertainty interval.
#'
#'   Selecting a topic estimate generates violin plots showing the p(s|k) distribution, split based
#'   on chosen binary covariate effects. The slider
#'   allows the user to threshold the number of points shown, based on their values in p(s|k).
#'   Highlighting points in the violin plots generates bar plots that show their abundances (or
#'   relative abundances) in the raw abundance table.
#'   }
#'
#' \subsection{Continuous}{
#'   Integrates the samples over topics p(t|s) and the topics
#'   over taxa p(k|t) distributions from the STM, binary and continuous covariate effects from the p(s|k) component, and
#'   their relationship with the raw taxonomic abundances. The covariate effects
#'   for each topic are shown as a scatterplot of posterior weights with error bars corresponding the
#'   global approximation of uncertainty. If the covariate chosen is binary, then the posterior regression
#'   weights with uncertainty intervals are shown. This is analogous to the mean difference between factor
#'   levels in the posterior predictive distribution. For continuous covariates, the points
#'   again represent the mean regression weights (i.e., the posterior slope estimate of the
#'   covariate). If, however, a spline or polynomial expansion was used, then the figure shows the posterior estimates of
#'   the standard deviation of the predicted topic probabilities from the posterior predictive distribution.
#'
#'   Selecting a topic estimate generates three panels. The top panel shows the posterior estimates of the
#'   selected continuous covariate. If binary covariates were present in the model formula, then the continuous effect
#'   given the binary covariate is shown as two regression lines, along with their corresponding uncertainty intervals.
#'   The points show the true p(k|s) values determined by the STM as a function of the selected continuous covariate.
#'   The middle panel then shows the raw abundances (or relative abundances) of most relavent taxa. Relavence can be
#'   control by adjusting \eqn{\lambda} where \deqn{r = \lambda x log p(t|k) + \lambda x log p(t|k)/p(x)}
#'   If binary
#'   covariates were provided in the model formula, selected split will split the regressions based on the selected
#'   covariate. Each figure overlays a linear best fit (red) and loess fit (red) to facilitate interpretation. The bottom
#'   panel shows these taxa combined.
#'   }
#'
#' \subsection{Functions}{
#'   Integrates the taxa over topics p(t|k) and gene functions over topics p(g|k) distributions,
#'   along with and the covariate effects from the p(s|k) component. The covariate effects
#'   for each topic are shown as a scatterplot of posterior weights with error bars corresponding the
#'   global approximation of uncertainty. If the covariate chosen is binary, then the posterior regression
#'   weights with uncertainty intervals are shown. This is analogous to the mean difference between factor
#'   levels in the posterior predictive distribution. For continuous covariates, the points
#'   again represent the mean regression weights (i.e., the posterior slope estimate of the
#'   covariate). If, however, a spline or polynomial expansion was used, then the figure shows the posterior estimates of
#'   the standard deviation of the predicted topic probabilities from the posterior predictive distribution.
#'   Colors indicate whether a given point was positive (red) or negative
#'   (blue) and did not enclose 0 at a user defined uncertainty interval.
#'
#'   The upper heatmap shows p(t|k), clustered via Wards method on a user chosen distance metric. Topics are ranked
#'   to right based on the weights from the aforementioned scatterplot. The lower heatmap shows the weights for the
#'   pathway-topic interaction from the multilevel Bayesian model. Positive and negative weight estimates that do
#'   not enclose zero at a chosen uncertainty level are marked with red and blue crosses, respectively. The pathway
#'   ordering is done via Wards method on Euclidean distance.
#'   Upon selected a cell within the pathway-topic heatmap, a table of genes is returned, ranking the genes in terms of
#'   abundance that belong to a given pathway-topic combination.
#'   }
#'
#' @references
#' Roberts, M.E., Stewart, B.M., Tingley, D., Lucas, C., Leder-Luis,
#' J., Gadarian, S.K., Albertson, B., & Rand, D.G. (2014). Structural topic
#' models for open-ended survey responses. Am. J. Pol. Sci. 58, 1064–1082.
#'
#' Sievert, C., & Shirley, K. (2014). LDAvis: A method for visualizing and
#' interpreting topics. Proc. Work. Interact. Lang. Learn. Vis. Interfaces.
#'
#' Zhao, T., & Liu., H. (2012) The huge Package for High-dimensional Undirected
#' Graph Estimation in R. Journal of Machine Learning Research.
#'
#' @seealso \code{\link[networkD3]{igraph_to_networkD3}}, \code{\link[huge]{huge}}, \code{\link[stm]{topicCorr}},
#' \code{\link[Rtsne]{Rtsne}}
#'
#' @examples
#' formula <- ~DIAGNOSIS
#' refs <- 'Not IBD'
#'
#' dat <- prepare_data(otu_table=GEVERS$OTU,rows_are_taxa=FALSE,tax_table=GEVERS$TAX,
#'                     metadata=GEVERS$META,formula=formula,refs=refs,
#'                     cn_normalize=TRUE,drop=TRUE)
#'
#' \dontrun{
#' vis(topic_effects,type='taxa')
#' vis(topic_effects,type='binary')
#' }
#'
#' formula <- ~PCDAI
#'
#' dat <- prepare_data(otu_table=GEVERS$OTU,rows_are_taxa=FALSE,tax_table=GEVERS$TAX,
#'                     metadata=GEVERS$META,formula=formula,refs=refs,
#'                     cn_normalize=TRUE,drop=TRUE)
#'
#' \dontrun{
#' vis(topic_effects,type='continuous')
#'
#' functions <- predict(topics,reference_path='/references/ko_13_5_precalculated.tab.gz')
#'
#' function_effects <- est(functions,level=3,
#'                         iters=500,method='hmc',
#'                         prior=c('laplace','t','laplace'))
#'
#' vis(function_effects,topic_effects)
#' }
#' @export
vis <- function(object,...) UseMethod('vis')

#' @rdname vis
#'
#' @param type Type of visualization to perform.
#' @param seed Seed for the random number generator to reproduce previous
#'   results.
#' @export
vis.effects <- function(object,topic_effects,type=c('taxa','binary','continuous','functions'),seed=object$seed$next_seed,...){

  set.seed(check_seed(seed))

  type <- match.arg(type)

  if (type == 'functions' | !missing(topic_effects)){

    objects <- list(object,topic_effects)
    types <- sapply(objects,function(x) attr(x,'type'))
    objects <- objects[match(types,c('functions','topics'))]

    if (any(sapply(objects,is.null))) stop('Must provide a topics effects-class and a functions effects-class.')

    object <- objects[[1]]
    object2 <- objects[[2]]
    class(object) <- 'functions'

    vis(object,object2,...)

  }else{

    class(object) <- type

    vis(object,...)

  }

}

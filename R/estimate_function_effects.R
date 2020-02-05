#' Estimate predicted function-topic effects
#'
#' Given within topic functional predictions, estimate the effects at a given
#' gene function category level. The effects correspond to a topic-gene category
#' interaction term after accounting for topic and gene category effects. The
#' model can be fit via either maximum likelihood or Hamiltonian MC.
#'
#' @param object (required) Ouput of \code{\link{predict.topics}}.
#' @param topics_subset Vector of topic indexes to be evaluated. Recommended to be < 25.
#' @param level Gene category level to evaluate. Defaults to 2.
#' @param method String indicating either ml or hmc. Defaults to hmc.
#' @param seed Seed for the random number generator to reproduce previous
#'   results.
#' @param verbose Logical flag to print progress information. Defaults to FALSE.
#' @param ... Additional arguments for methods.
#'
#' @return An object of class effects containing
#' \describe{
#' \item{model}{List containing the parameters, fit, and summary.}
#' \item{gene_table}{Dataframe containing the formatted predicted gene information
#' from \code{\link{predict.topics}}.}
#' }
#'
#' @details
#' The functional effects are estimated via a multilevel Bayesian negative binomial regression model.
#' Topic and pathway level effects are estimated, as well as topic-pathway interactions. The model
#' has the following form:
#' \deqn{\theta_{i} = \mu + \beta_{w} + \beta_{k} + \beta_{w,k}}
#' \deqn{y_{i} ~ NB(\theta_{i},\phi)}
#  gene counts y are distributed by a negative binomial distribution with mean \eqn{\theta},
#' where \eqn{\mu} is the intercept and each \eqn{\beta} term represents the weight for pathway level,
#' topic, and pathway level-topic interaction, respectively; \eqn{\phi} represents the dispersion
#' parameter.
#'
#' \subsection{HMC}{Hamiltonian MC is performed via Stan. By default, student-t priors with degrees of
#' freedom set at 7 are placed on all regression weights, with variance terms distributed by half normal
#' priors. The intercept \eqn{\mu} is given a normal prior with fixed variance. Lastly, \eqn{\phi} is
#' given an \eqn{exponential(.5)} prior. The priors placed on the regression weights can be changed by
#' the user to either normal, t-family, or laplace (double exponential) priors if a sparse solution is
#' desired. For the latter, each variance term is given an additional regularization parameter
#' \eqn{\lambda} which in turn is distributed by a \eqn{chi-squared(1)} distribution.
#'
#' Unless a set of initialization values are provided by the user, or the user chooses to select a random
#' initialization procedure, initial values are set at the maximum likelihood estimate via
#' \code{\link[lme4]{glmer.nb}}, but at a far smaller number of iterations than had the user chosen
#' ML as his or her estimation method.}
#'
#' \subsection{ML}{Maximum likelihood estimation is performed via \code{\link[lme4]{glmer.nb}}. For
#' deeper level functional categories, the model may fail to converge, even with a substantial
#' number of iterations. In such a case, the model estimates are returns so the user can perform
#' HMC, but by initializing at these ML values.}
#'
#' @seealso \code{\link[lme4]{glmer.nb}} \code{\link[rstan]{stan}} \code{\link{resume}}
#'
#' @references
#' Bates, D., Maechler, M., Bolker, B., and Walker, S. (2015).
#' Fitting Linear Mixed-Effects Models Using lme4. Journal of
#' Statistical Software, 67(1), 1-48. doi:10.18637/jss.v067.i01.
#'
#' Gelman, A. and Hill, J. (2006). Data Analysis Using Regression and
#' Multilevel/Hierarchical Models. Cambridge University Press; 1 edition.
#'
#' Stan Development Team. 2016. RStan: the R interface to Stan.
#' http://mc-stan.org
#'
#' Stan Development Team. 2016. Stan Modeling Language Users Guide and
#' Reference Manual, Version 2.14.0. http://mc-stan.org
#'
#' @examples
#' formula <- ~DIAGNOSIS
#' refs <- 'CD'
#'
#' dat <- prepare_data(otu_table=GEVERS$OTU,rows_are_taxa=FALSE,tax_table=GEVERS$TAX,
#'                     metadata=GEVERS$META,formula=formula,refs=refs,
#'                     cn_normalize=TRUE,drop=TRUE)
#'
#' \dontrun{
#' topics <- find_topics(dat,K=15)
#'
#' functions <- predict(topics,reference_path='/references/ko_13_5_precalculated.tab.gz')
#' function_effects <- est(functions,level=3,
#'                         iters=500,method='hmc',
#'                         prior=c('laplace','t','laplace'))
#' }
#'
#' @aliases est_functions est.functions
#'
#' @export

est.functions <- function(object,topics_subset,level=2,method=c('hmc','ml'),
                          seed=object$seeds$next_seed,verbose=FALSE,...){

  set.seed(check_seed(seed))
  mod_seed <- sample.int(.Machine$integer.max,1)
  next_seed <- sample.int(.Machine$integer.max,1)

  method <- match.arg(method)

  fxn_table <- object$fxn_table
  fxn_meta <- object$fxn_meta

  if (missing(topics_subset) & nrow(fxn_table) <= 25) topics_subset <- 1:nrow(fxn_table)

  fxn_table <- fxn_table[topics_subset,]
  fxn_table <- fxn_table[,colSums(fxn_table) > 0]

  fxn_meta <- lapply(fxn_meta,function(x) x[colnames(fxn_table)])

  object$fxn_table <- fxn_table
  object$fxn_meta <- fxn_meta

  gene_table <- list(gene_table=format_gene_table(object,level=level))

  class(gene_table) <- method
  mm <- est(gene_table,seed=mod_seed,verbose=verbose,...)

  out <- list(model=mm,gene_table=gene_table$gene_table,seeds=list(seed=seed,model_seed=mod_seed,next_seed=next_seed))
  class(out) <- 'effects'
  attr(out,'type') <- 'functions'
  attr(out,'method') <- method

  return(out)

}

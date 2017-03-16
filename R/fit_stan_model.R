#' @importFrom lme4 glmer.nb glmer
#' @importFrom rstan stan
#' @importFrom stats4 summary
NULL

#' Estimate function effects via HMC
#'
#' Given within topic functional predictions, estimate the effects at a given
#' gene function category level via HMC. The effects correspond to a topic-gene
#' category interaction term after accounting for topic and gene category
#' effects.
#'
#' @param gene_table A gene_table formatted via \code{\link{format_gene_table}}.
#' @param inits List of values for parameter initialization. If omitted, values
#'   are generated via \code{\link{glmer.nb}}
#' @param iters Number of iterations for HMC. Defaults to 1000.
#' @param chains (optional) Number of chains for HMC. Defaults to 1.
#' @param return_fit (optional) Logical to return the rstan data and fit.
#' @param verbose (optional) Logical to print progress.
#'
#' @return A list containing
#'
#' \item{summary}{Rstan output that includes the coefficient mean estimates, uncertainty
#' intervals, standard deviations, effective sample size, Rhat statistics}
#' \item{fit}{Rstan fit (if return_fit=TRUE)}
#' \item{data}{Rstan input data (if return_fit=TRUE)}
#' \item{flagged}{Parameter estimates with Rhat statistics > 1.1}
#' \item{inits}{If parameter estimates are flagged, a list of inits to initialize Rstan
#' with more iterations.}
#' @export

fit_stan_model <- function(gene_table,inits,iters=1000,chains=1,return_fit=FALSE,verbose=FALSE){

  stan_dat <- list(N=nrow(gene_table),
                   J=length(unique(gene_table$pw)),
                   K=length(unique(gene_table$topic)),
                   I=length(levels(as.factor(gene_table$pw):as.factor(gene_table$topic))),
                   pw_full=gene_table$pw,
                   pw=as.integer(as.factor(gene_table$pw)),
                   topic_full=gene_table$topic,
                   topic=gene_table$topic,
                   pwxtopic_full=as.factor(gene_table$pw):as.factor(gene_table$topic),
                   pwxtopic=as.integer(as.factor(gene_table$pw):as.factor(gene_table$topic)),
                   y=gene_table$count)

  stan_code <- '
  data{
    int<lower=1> N;
    int<lower=1> J;
    int<lower=1> K;
    int<lower=1> I;

    int<lower=0> y[N];

    int<lower=1,upper=K> topic[N];
    int<lower=1,upper=J> pw[N];
    int<lower=1,upper=I> pwxtopic[N];
  }
  parameters{
    vector[J] b_pw;
    vector[K] b_topic;
    vector[I] b_pwxtopic;

    real mu;

    real<lower=0,upper=100> b_pw_sigma;
    real<lower=0,upper=100> b_topic_sigma;
    real<lower=0,upper=100> b_pwxtopic_sigma;

    real<lower=0> phi;
  }
  model{
    real theta;

    phi ~ exponential(1.25);

    mu ~ normal(0,5);
    b_pw ~ normal(0,b_pw_sigma);
    b_topic ~ normal(0,b_topic_sigma);
    b_pwxtopic ~ normal(0,b_pwxtopic_sigma);

    b_pw_sigma ~ normal(0,2.5);
    b_topic_sigma ~ normal(0,2.5);
    b_pwxtopic_sigma ~ normal(0,2.5);

    for (n in 1:N){
      theta = mu + b_pw[pw[n]] + b_topic[topic[n]] + b_pwxtopic[pwxtopic[n]];
      y[n] ~ neg_binomial_2_log(theta,phi);
    }
  }
  generated quantities{
    real theta;
    vector[N] yhat;
    vector[N] log_lik;

    for (n in 1:N){
      theta = mu + b_pw[pw[n]] + b_topic[topic[n]] + b_pwxtopic[pwxtopic[n]];
      yhat[n] = neg_binomial_2_log_rng(theta,phi);
      log_lik[n] = neg_binomial_2_log_lpmf(y[n]|theta,phi);
    }
  }
  '

  stan_table <- data.frame(y=stan_dat$y,
                           pw=stan_dat$pw,
                           topic=stan_dat$topic,
                           pwxtopic=stan_dat$pwxtopic)



  if (missing(inits)){

    if (verbose) cat('Generating initial values.\n')

    mm_init <- glmer.nb(y ~ (1|pw) + (1|topic) + (1|pwxtopic),
                        data=stan_table,
                        verbose=verbose)

    inits <- lapply(seq_len(chains),function(x) list(mu=fixef(mm_init),
                                                     phi=getME(mm_init,'glmer.nb.theta'),
                                                     b_pw=unlist(ranef(mm_init)$pw),
                                                     b_topic=unlist(ranef(mm_init)$topic),
                                                     b_pwxtopic=unlist(ranef(mm_init)$pwxtopic)))

  }else if (class(inits) == 'glmerMod'){

    inits <- lapply(seq_len(chains),function(x) list(mu=fixef(inits),
                                                     phi=getME(inits,'glmer.nb.theta'),
                                                     b_pw=unlist(ranef(inits)$pw),
                                                     b_topic=unlist(ranef(inits)$topic),
                                                     b_pwxtopic=unlist(ranef(inits)$pwxtopic)))

  }else if (class(inits) == 'list'){

    inits <- inits[names(inits) %in% c('mu','phi','b_pw','b_topic','b_pwxtopic')]

    inits <- lapply(seq_len(chains),function(x) inits)

  }else if (is.null(inits) | tolower(inits) == 'none' | is.na(inits) | tolower(inits) == 'random'){

    inits <- 'random'

  }else{

    stop('inits must be a list, glmer.nb object (glmerMod), or NULL.')

  }



  if (chains > 1){
    if (verbose) cat('Preparing parallelization.\n')
    options_old <- options()

    on.exit(options(options_old),add=TRUE)

    rstan::rstan_options(auto_write=TRUE)
    options(mc.cores=chains)
  }



  if (verbose) cat('Fitting model via HMC.\n')

  fit <- stan(model_code=stan_code,data=stan_dat,
                     pars=c('theta'),include=FALSE,
                     init=inits,
                     iter=iters,chains=chains,
                     verbose=verbose)

  if (verbose) cat('Extracting summary.\n')

  summary_pars <- c('mu','phi','b_pw_sigma','b_topic_sigma','b_pwxtopic_sigma','b_pw','b_topic','b_pwxtopic','yhat')
  extract_summary <- vector(mode='list',length=length(summary_pars))
  names(extract_summary) <- summary_pars
  for (i in seq_along(extract_summary)){

    if (grepl('^b\\_',summary_pars[i])){

      extract_summary_tmp <- summary(fit,pars=summary_pars[i])[['summary']]
      par_name <- gsub('^b\\_','',summary_pars[i])
      lookup_table <- unique(cbind(as.character(stan_dat[[par_name]]),as.character(stan_dat[[paste0(par_name,'_full')]])))
      rownames(extract_summary_tmp) <- lookup_table[,2][order(as.integer(lookup_table[,1],decreasing=TRUE))]

      if (par_name == 'pwxtopic'){
        par_name_tmp <- do.call('rbind',strsplit(rownames(extract_summary_tmp),'\\:'))
        colnames(par_name_tmp) <- unlist(strsplit(par_name,'x'))
        extract_summary_tmp <- cbind(par_name_tmp,as.data.frame(extract_summary_tmp))
      }else{
        par_name_tmp <- matrix(rownames(extract_summary_tmp),ncol=1)
        colnames(par_name_tmp) <- par_name
        extract_summary_tmp <- cbind(par_name_tmp,as.data.frame(extract_summary_tmp))
      }

      extract_summary[[i]] <- extract_summary_tmp

    }else{

      extract_summary[[i]] <- summary(fit,pars=summary_pars[i])[['summary']]

    }

  }

  out <- list(summary=extract_summary)

  if (return_fit){
    out[['fit']] <- fit
    out[['data']] <- stan_dat
  }

  rhat <- summary(fit,pars=c('mu','phi','b_pw','b_topic','b_pwxtopic'))[['summary']][,'Rhat'] > 1.1
  rhat_count <- sum(rhat)
  if (rhat_count > 0){

    warning(sprintf('%s parameters with Rhat > 1.1. Consider more iterations.',rhat_count))

    out[['flagged']] <- names(which(rhat))
    out[['inits']] <- list(mu=summary(fit,pars='mu')[['summary']][,'mean'],
                           phi=summary(fit,pars='phi')[['summary']][,'mean'],
                           b_pw=summary(fit,pars='b_pw')[['summary']][,'mean'],
                           b_topic=summary(fit,pars='b_topic')[['summary']][,'mean'],
                           b_pwxtopic=summary(fit,pars='b_pwxtopic')[['summary']][,'mean'])

  }

  out[['sampler']] <- rstan::get_sampler_params(fit)

  return(out)

}

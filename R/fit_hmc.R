#' @rdname est.functions
#'
#' @param inits List of values for parameter initialization. If omitted, values
#'   are generated via \code{\link{glmer.nb}}
#' @param prior Prior to be placed on covariate weights. Choices include student-t,
#' normal, and laplace. Defaults to student-t.
#' @param t_df Degrees of freedom for student-t priors. Defaults to 7.
#' @param iters Number of iterations for for fitting. Defaults to 300 and 100 for
#' HMC and ML, respectively.
#' @param warmup For HMC, proportion of iterations devoted to warmup. Defaults to
#' iters/2.
#' @param chains For HMC, number of independent chains. Defaults to 1.
#' @param cores For HMC, number of cores to parallelize chains. Defaults to 1.
#' @param return_summary Logical flag to return results summary. Defaults to TRUE.
#'
#' @export

est.hmc <- function(object,inits,prior=c('t','normal','laplace'),t_df=c(7,7,7),iters=300,warmup=iters/2,
                    chains=1,cores=1,seed=sample.int(.Machine$integer.max,1),
                    return_summary=TRUE,verbose=FALSE,...){

  set.seed(check_seed(seed))
  mod_seed <- sample.int(.Machine$integer.max,1)
  next_seed <- sample.int(.Machine$integer.max,1)

  gene_table <- object$gene_table

  prior <- match.arg(prior,several.ok=TRUE)

  lambda_params <- c('','','')
  lambda_rhat <- vector(mode='character')
  priors <- c('b_pw ~ student_t(7,0,b_pw_sigma);\n',
              'b_topic ~ student_t(7,0,b_topic_sigma);\n',
              'b_pwxtopic ~ student_t(7,0,b_pwxtopic_sigma);\n')
  priors_level <- c('pw','topic','pwxtopic')
  for (i in seq_along(prior)){
    if (prior[i] == 't'){
      priors[i] <- sprintf('b_%s ~ student_t(%s,0,b_%s_sigma);\n',priors_level[i],t_df[1],priors_level[i])
      t_df <- t_df[-1]
    }else if (prior[i] == 'normal'){
      priors[i] <- sprintf('b_%s ~ normal(0,b_%s_sigma);\n',priors_level[i],priors_level[i])
    }else if (prior[i] == 'laplace'){
      lambda_params[i] <- sprintf('real lambda_%s;\n',priors_level[i],priors_level[i])
      lambda_rhat <- c(lambda_rhat,sprintf('lambda_%s',priors_level[i]))
      priors[i] <- sprintf('lambda_%s ~ chi_square(1);\nb_%s ~ double_exponential(0,b_%s_sigma/lambda_%s);\n',
                           priors_level[i],priors_level[i],priors_level[i],priors_level[i])
    }
  }

  if (verbose) cat(sprintf('Setting regression coefficient priors to\n%s\n',paste0(priors,collapse='')))

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

  stan_code <- sprintf('
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

    %s // lambda_pw
    %s // lambda_topic
    %s // lambda_pwxtopic

    real<lower=0> phi;
  }
  model{
    vector[N] theta;

    phi ~ exponential(1.25);

    mu ~ normal(0,5);

    %s // b_pw
    %s // b_topic
    %s // b_pwxtopic

    b_pw_sigma ~ normal(0,2.5);
    b_topic_sigma ~ normal(0,2.5);
    b_pwxtopic_sigma ~ normal(0,2.5);

    for (n in 1:N)
      theta[n] = mu + b_pw[pw[n]] + b_topic[topic[n]] + b_pwxtopic[pwxtopic[n]];

    y ~ neg_binomial_2_log(theta,phi);

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
    ',lambda_params[1],lambda_params[2],lambda_params[3],priors[1],priors[2],priors[3])

  stan_table <- data.frame(y=stan_dat$y,
                           pw=stan_dat$pw,
                           topic=stan_dat$topic,
                           pwxtopic=stan_dat$pwxtopic,
                           stringsAsFactors=TRUE)

  if (missing(inits)){

    if (verbose) cat('Generating initial values via ML.\n')

    mm_init <- suppressWarnings(
      glmer.nb(y ~ (1|pw) + (1|topic) + (1|pwxtopic),
                        data=stan_table,
                        verbose=verbose,
                        control=glmerControl(calc.derivs=FALSE,
                                             optCtrl=list(maxfun=50))) # check if ok for level 3
    )

    inits <- lapply(seq_len(chains),function(x) list(mu=fixef(mm_init),
                                                     phi=getME(mm_init,'glmer.nb.theta'),
                                                     b_pw=unlist(ranef(mm_init)$pw),
                                                     b_topic=unlist(ranef(mm_init)$topic),
                                                     b_pwxtopic=unlist(ranef(mm_init)$pwxtopic)))

  }else if (class(inits)[1] == 'glmerMod'){

    inits <- lapply(seq_len(chains),function(x) list(mu=fixef(inits),
                                                     phi=getME(inits,'glmer.nb.theta'),
                                                     b_pw=unlist(ranef(inits)$pw),
                                                     b_topic=unlist(ranef(inits)$topic),
                                                     b_pwxtopic=unlist(ranef(inits)$pwxtopic)))

  }else if (class(inits)[1] == 'list'){

    inits <- inits[names(inits) %in% c('mu','phi','b_pw','b_topic','b_pwxtopic')]

    inits <- lapply(seq_len(chains),function(x) inits)

  }else if (is.null(inits) | tolower(inits) == 'none' | is.na(inits) | tolower(inits) == 'random'){

    inits <- 'random'

  }else{

    stop('inits must be a list, glmer.nb object (glmerMod), or NULL.')

  }

  if (cores > 1){
    if (verbose) cat('Preparing parallelization.\n')
    options_old <- options()

    on.exit(options(options_old),add=TRUE)

    rstan::rstan_options(auto_write=TRUE)
    options(mc.cores=cores)
  }

  if (verbose) cat('Fitting model via HMC.\n')

  fit <- rstan::stan(model_code=stan_code,data=stan_dat,
                     pars=c('theta'),include=FALSE,
                     init=inits,
                     warmup=warmup,
                     iter=iters,
                     chains=chains,cores=cores,
                     seed=mod_seed,
                     verbose=verbose)

  out <- list()
  summary_pars <- c('mu','phi',lambda_rhat,'b_pw_sigma','b_topic_sigma','b_pwxtopic_sigma','b_pw','b_topic','b_pwxtopic','yhat')
  out[['pars']] <- summary_pars
  out[['fit']] <- fit
  out[['data']] <- stan_dat
  out[['inits']] <- list(orig=inits)
  out[['sampler']] <- rstan::get_sampler_params(fit)


  if (return_summary){

    if (verbose) cat('Extracting summary (this often takes some time).\n')
    out[['summary']] <- extract_stan_summary(fit,stan_dat,summary_pars)
    rhat_pars <- c('mu','phi',lambda_rhat,'b_pw','b_topic','b_pwxtopic','b_pw_sigma','b_topic_sigma','b_pwxtopic_sigma')
    rhat <- stats4::summary(fit,pars=rhat_pars)[['summary']][,'Rhat'] > 1.1
    rhat_count <- sum(rhat,na.rm=TRUE)
    if (rhat_count > 0){
      warning(sprintf('%s parameters with Rhat > 1.1. Consider more iterations.',rhat_count))
      out[['flagged']] <- names(which(rhat))
    }
  }

  out[['seeds']] <- list(seed=seed,mod_seed=mod_seed,next_seed=next_seed)

  return(out)

}

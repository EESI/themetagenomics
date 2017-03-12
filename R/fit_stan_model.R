fit_stan_model <- function(gene_table,iters=1000,chains=1){

  stan_dat <- list(N=nrow(gene_table),
                   J=length(unique(gene_table$pw)),
                   K=length(unique(gene_table$topic)),
                   I=length(levels(as.factor(gene_table$pw):as.factor(gene_table$topic))),
                   pw_full=gene_table$pw,
                   pw=as.integer(as.factor(gene_table$pw)),
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

  if (verbose) cat('Generating initial values.\n')

  mm_init <- lme4::glmer.nb(y ~ (1|pw) + (1|topic) + (1|pwxtopic),
                            data=stan_table,
                            verbose=verbose,
                            ...)

  inits <- lapply(seq_len(chains),function(x) list(mu=lme4::fixef(mm_init),
                                                   phi=lme4::getME(mm_init,"glmer.nb.theta"),
                                                   b_pw=unlist(lme4::ranef(mm_init)$pw),
                                                   b_topic=unlist(lme4::ranef(mm_init)$topic),
                                                   b_pwxtopic=unlist(lme4::ranef(mm_init)$pwxtopic)))

  if (chains > 1){
    if (verbose) cat('Preparing parallelization.\n')
    options_old <- options()
    rstan_options_old <- rstan_options()

    on.exit(options(options_old),add=TRUE)
    on.exit(rstan_options(rstan_options_old),add=TRUE)

    rstan_options(auto_write=TRUE)
    options(mc.cores = chains)
  }

  if (verbose) cat('Fitting model via Stan.\n')

  fit <- rstan::stan(model_code=stan_code,data=stan_dat,
                     pars=c('theta'),include=FALSE,
                     init=inits,
                     iter=iters,chains=chains,
                     verbose=verbose)

  return(list(fit=fit,data=stan_dat))

}

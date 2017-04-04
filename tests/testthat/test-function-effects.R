context("test function effects")

Sys.setenv("R_TESTS" = "")

test_that('picrust results gives similar results for ml and hmc',{

  suppressWarnings({
    ml <- est(GEVERS$functions,level=2,iters=100,method='ml')
    hmc <- est(GEVERS$functions,level=2,iters=50,chains=1,
               return_summary=TRUE,
               prior=c('t','t','normal'),
               t_df=c(7,7),seed=123)
    hmc_summary <- extract(hmc)
    resume_orig <- resume(hmc,init_type='orig',iters=50,chains=1,return_summary=TRUE,seed=123)
    resume_last <- resume(hmc,init_type='last',iters=50,chains=1,return_summary=TRUE,seed=123)
  })

  overlap <- c('mu','phi','b_pw_sigma','b_topic_sigma','b_pwxtopic_sigma','b_pw','b_topic','b_pwxtopic','yhat')

  expect_warning(extract(hmc))
  expect_identical(unlist(lapply(hmc$model$summary[overlap],rownames)),
                   unlist(lapply(ml$model$summary[overlap],rownames)))
  expect_identical(unlist(lapply(hmc$model$summary[overlap],rownames)),
                   unlist(lapply(hmc_summary$summary[overlap],rownames)))
  expect_identical(unlist(lapply(hmc$model$summary[overlap],rownames)),
                   unlist(lapply(resume_orig$model$summary[overlap],rownames)))
  expect_identical(unlist(lapply(hmc$model$summary[overlap],rownames)),
                   unlist(lapply(resume_last$model$summary[overlap],rownames)))

  suppressWarnings({
    ml <- est(GEVERS$functions,level=3,iters=5,method='ml')
    hmc <- est(GEVERS$functions,level=3,iters=50,chains=1,
               return_summary=TRUE,
               prior=c('t','t','normal'),
               t_df=c(7,7),seed=123)
  })

  overlap <- c('mu','phi','b_pw_sigma','b_topic_sigma','b_pwxtopic_sigma','b_pw','b_topic','b_pwxtopic','yhat')

  expect_identical(unlist(lapply(hmc$model$summary[overlap],rownames)),
                   unlist(lapply(ml$model$summary[overlap],rownames)))

})

test_that('tax4fun results gives similar results for ml and hmc',{

  suppressWarnings({
    ml <- est(DAVID$functions,level=2,iters=100,method='ml')
    hmc <- est(DAVID$functions,level=2,iters=50,chains=1,
               return_summary=TRUE,
               prior=c('t','t','normal'),
               t_df=c(7,7),seed=123)
    hmc_summary <- extract(hmc)
    resume_orig <- resume(hmc,init_type='orig',iters=50,chains=1,return_summary=TRUE,seed=123)
    resume_last <- resume(hmc,init_type='last',iters=50,chains=1,return_summary=TRUE,seed=123)
  })

  overlap <- c('mu','phi','b_pw_sigma','b_topic_sigma','b_pwxtopic_sigma','b_pw','b_topic','b_pwxtopic','yhat')

  expect_warning(extract(hmc))
  expect_identical(unlist(lapply(hmc$model$summary[overlap],rownames)),
                   unlist(lapply(ml$model$summary[overlap],rownames)))
  expect_identical(unlist(lapply(hmc$model$summary[overlap],rownames)),
                   unlist(lapply(hmc_summary$summary[overlap],rownames)))
  expect_identical(unlist(lapply(hmc$model$summary[overlap],rownames)),
                   unlist(lapply(resume_orig$model$summary[overlap],rownames)))
  expect_identical(unlist(lapply(hmc$model$summary[overlap],rownames)),
                   unlist(lapply(resume_last$model$summary[overlap],rownames)))

  suppressWarnings({
    ml <- est(DAVID$functions,level=3,iters=5,method='ml')
    hmc <- est(DAVID$functions,level=3,iters=50,chains=1,
               return_summary=TRUE,
               prior=c('t','t','normal'),
               t_df=c(7,7),seed=123)
  })

  overlap <- c('mu','phi','b_pw_sigma','b_topic_sigma','b_pwxtopic_sigma','b_pw','b_topic','b_pwxtopic','yhat')

  expect_identical(unlist(lapply(hmc$model$summary[overlap],rownames)),
                   unlist(lapply(ml$model$summary[overlap],rownames)))

})

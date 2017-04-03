context("test function effects")

test_that('picrust results gives similar results for ml and hmc',{

  # expect_warning(est(GEVERS$functions,level=2,iters=5,method='ml'))
  expect_warning(est(GEVERS$functions,level=2,iters=50,chains=1,
                     return_summary=TRUE,
                     prior=c('t','t','normal'),
                     t_df=c(7,7),seed=123))

})

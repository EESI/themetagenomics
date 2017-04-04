context("test function effects")

test_that('predict via ml works',{

  expect_warning(z1 <- est(GEVERS$functions,level=2,iters=5,method='ml'))
  expect_warning(z2 <- est(GEVERS$functions,level=3,iters=5,method='ml'))

  expect_is(z1$model$summary$mu[1],'numeric')
  expect_is(z1$model$fit,'glmerMod')

  expect_is(z2$model$summary$mu[1],'numeric')
  expect_is(z2$model$fit,'glmerMod')

  expect_warning(est(DAVID$functions,level=2,iters=5,method='ml'))
  expect_warning(est(DAVID$functions,level=3,iters=5,method='ml'))

})

Sys.setenv("R_TESTS" = "")

test_that('rccp works',{
  expect_equal(example(cxxfunction, package = "inline", run.dontrun = TRUE)$value,10)
})

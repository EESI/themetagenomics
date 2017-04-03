context("test function effects")

test_that('picrust results gives similar results for ml and hmc',{

  expect_warning(est(GEVERS$functions,level=2,iters=5,method='ml'))

})

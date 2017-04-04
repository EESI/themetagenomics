context("test function effects")

Sys.setenv("R_TESTS" = "")

test_that('rccp works',{
  expect_equal(example(cxxfunction, package = "inline", run.dontrun = TRUE)$value,10)
})

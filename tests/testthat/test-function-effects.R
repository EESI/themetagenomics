context("test function effects")

Sys.setenv("R_TESTS" = "")

test_that('rccp works',{
  expect_equal(example(cxxfunction, package = "inline", run.dontrun = TRUE)$value,10)
})

test_that('glmer works',{
  DAT <- readRDS(system.file('testdata','otufuncdata.rds',package='themetagenomics'))
  ml <- suppressWarnings(est(DAT,level=2,iters=5,method='ml'))
  expect_is(ml,'effects')
})

test_that('stan works',{
  DAT <- readRDS(system.file('testdata','otufuncdata.rds',package='themetagenomics'))
  hmc <- suppressMessages(suppressWarnings(est(DAT,level=2,iters=50,chains=1,method='hmc')))
  expect_is(hmc,'effects')
})

context("function effects")

Sys.setenv("R_TESTS" = "")

test_that('picrust results gives similar results for ml and hmc',{

  skip_on_cran()

  DAT <- readRDS(system.file('testdata','otufuncdata.rds',package='themetagenomics'))
  overlap <- c('mu','phi','b_pw_sigma','b_topic_sigma','b_pwxtopic_sigma','b_pw','b_topic','b_pwxtopic','yhat')

  ml <- suppressMessages(suppressWarnings(est(DAT,level=2,iters=5,method='ml')))
  hmc <- suppressMessages(suppressWarnings(est(DAT,level=2,iters=50,chains=1,
             return_summary=TRUE,
             prior=c('t','t','normal'),
             t_df=c(7,7),seed=123)))
  hmc_summary <- suppressMessages(suppressWarnings(extract(hmc)))
  resume_orig <- suppressMessages(suppressWarnings(resume(hmc,init_type='orig',iters=50,chains=1,return_summary=TRUE)))
  resume_last <- suppressMessages(suppressWarnings(resume(hmc,init_type='last',iters=50,chains=1,return_summary=TRUE)))

  expect_identical(unlist(lapply(hmc$model$summary[overlap],rownames)),
                   unlist(lapply(ml$model$summary[overlap],rownames)))
  expect_identical(unlist(lapply(hmc$model$summary[overlap],rownames)),
                   unlist(lapply(hmc_summary$summary[overlap],rownames)))
  expect_identical(unlist(lapply(hmc$model$summary[overlap],rownames)),
                   unlist(lapply(resume_orig$model$summary[overlap],rownames)))
  expect_identical(unlist(lapply(hmc$model$summary[overlap],rownames)),
                   unlist(lapply(resume_last$model$summary[overlap],rownames)))

})

test_that('format_gene_table works for level 3',{

  DAT <- readRDS(system.file('testdata','otufuncdata.rds',package='themetagenomics'))
  gene_table2 <- format_gene_table(DAT,level=2)
  gene_table3 <- format_gene_table(DAT,level=3)

  expect_is(gene_table2,'data.frame')
  expect_true(max(gene_table2$count) > 0)
  expect_true(min(gene_table2$count) == 0)

  expect_is(gene_table3,'data.frame')
  expect_true(max(gene_table3$count) > 0)
  expect_true(min(gene_table3$count) == 0)

  expect_true(length(unique(gene_table3$pw)) > length(unique(gene_table2$pw)))

  DAT <- readRDS(system.file('testdata','seqfuncdata.rds',package='themetagenomics'))
  gene_table2 <- format_gene_table(DAT,level=2)
  gene_table3 <- format_gene_table(DAT,level=3)

  expect_is(gene_table2,'data.frame')
  expect_true(max(gene_table2$count) > 0)
  expect_true(min(gene_table2$count) == 0)

  expect_is(gene_table3,'data.frame')
  expect_true(max(gene_table3$count) > 0)
  expect_true(min(gene_table3$count) == 0)

  expect_true(length(unique(gene_table3$pw)) > length(unique(gene_table2$pw)))

})

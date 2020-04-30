context('estimate topic effects')

test_that('est.topics returns correct results for different formulae',{

  skip_on_cran()
  skip_on_travis()

  DAT <- readRDS(system.file('testdata','otudata.rds',package='themetagenomics'))

  set.seed(23)

  x <- prepare_data(otu_table=DAT$OTU,
                    rows_are_taxa=FALSE,
                    tax_table=DAT$TAX,
                    metadata=DAT$META,
                    formula=~DIAGNOSIS,
                    refs='Not IBD',
                    cn_normalize=TRUE,
                    drop=TRUE)

  y <- find_topics(x,K=5,init_type='Spectral',tol=1e-03)

  z1 <- est(y)

  expect_identical(colnames(est(y)[[1]][[1]][[1]]),c('estimate','10%','90%'))
  expect_identical(colnames(est(y,ui_level=.95)[[1]][[1]][[1]]),c('estimate','2.5%','97.5%'))
  expect_identical(y$modelframe,z1$modelframe)

  expect_error(est(y,metadata=DAT$META[1:7,],formula=~PCDAI))

  x <- prepare_data(otu_table=DAT$OTU,
                    rows_are_taxa=FALSE,
                    tax_table=DAT$TAX,
                    metadata=DAT$META,
                    formula=~PCDAI,
                    refs='Not IBD',
                    cn_normalize=TRUE,
                    drop=TRUE)

  y <- find_topics(x,K=5,init_type='Spectral',tol=1e-03)
  expect_identical(y$modelframe,est(y)$modelframe)

  expect_error(est(y,metadata=DAT$META[1:7,],formula=~DIAGNOSIS))
  expect_warning(est(y,metadata=DAT$META,formula=~DIAGNOSIS))

  z2 <- est(y,metadata=DAT$META,formula=~DIAGNOSIS,refs='Not IBD')

  expect_identical(z1$modelframe[rownames(DAT$META)[!is.na(DAT$META$PCDAI)],],
                   z2$modelframe[rownames(DAT$META)[!is.na(DAT$META$PCDAI)],])

  DAT <- readRDS(system.file('testdata','seqdata.rds',package='themetagenomics'))

  set.seed(423)

  x <- prepare_data(otu_table=DAT$ABUND,
                    rows_are_taxa=FALSE,
                    tax_table=DAT$TAX,
                    metadata=DAT$META,
                    formula=~Site + Day,
                    refs='UBERON:saliva',
                    cn_normalize=FALSE,
                    drop=TRUE)

  y <- find_topics(x,K=5,init_type='Spectral')
  z1 <- est(y)

  expect_identical(colnames(est(y)[[1]][[1]][[1]]),c('estimate','10%','90%'))
  expect_identical(colnames(est(y,ui_level=.95)[[1]][[1]][[1]]),c('estimate','2.5%','97.5%'))
  expect_identical(y$modelframe,z1$modelframe)

  x <- prepare_data(otu_table=DAT$ABUND,
                    rows_are_taxa=FALSE,
                    tax_table=DAT$TAX,
                    metadata=DAT$META,
                    formula=~Day + Site,
                    refs='UBERON:saliva',
                    cn_normalize=FALSE,
                    drop=TRUE)

  y <- find_topics(x,K=5,init_type='Spectral',tol=1e-03)
  z2 <- est(y)

  expect_true(mean(abs(z1$topic_effects$`SiteUBERON:feces`$est[,1]-z2$topic_effects$`SiteUBERON:feces`$est[,1])) < .1)
  expect_true(mean(abs(z1$topic_effects$Day$est[,1]-z2$topic_effects$Day$est[,1])) < .1)

  set.seed(23)

  x <- prepare_data(otu_table=DAT$ABUND,
                    rows_are_taxa=FALSE,
                    tax_table=DAT$TAX,
                    metadata=DAT$META,
                    formula=~Site + s(Day),
                    refs='UBERON:saliva',
                    cn_normalize=FALSE,
                    drop=TRUE)

  y <- find_topics(x,K=5,init_type='Spectral',tol=1e-03)
  z1 <- est(y)

  expect_identical(colnames(est(y)[[1]][[1]][[1]]),c('estimate','10%','90%'))
  expect_identical(colnames(est(y,ui_level=.95)[[1]][[1]][[1]]),c('estimate','2.5%','97.5%'))
  expect_identical(y$modelframe,z1$modelframe)

  x <- prepare_data(otu_table=DAT$ABUND,
                    rows_are_taxa=FALSE,
                    tax_table=DAT$TAX,
                    metadata=DAT$META,
                    formula=~s(Day) + Site,
                    refs='UBERON:saliva',
                    cn_normalize=FALSE,
                    drop=TRUE)

  y <- find_topics(x,K=5,init_type='Spectral',tol=1e-03)
  z2 <- est(y)

  expect_true(mean(abs(z1$topic_effects$`SiteUBERON:feces`$est[,1]-z2$topic_effects$`SiteUBERON:feces`$est[,1])) < .1)
  expect_true(mean(abs(z1$topic_effects$Day$est[,1]-z2$topic_effects$Day$est[,1])) < .1)

  set.seed(23)

  x <- prepare_data(otu_table=DAT$ABUND,
                    rows_are_taxa=FALSE,
                    tax_table=DAT$TAX,
                    metadata=DAT$META,
                    formula=~Multi + s(Day),
                    refs='1',
                    cn_normalize=FALSE,
                    drop=TRUE)

  y <- find_topics(x,K=5,init_type='Spectral',tol=1e-03)
  z1 <- est(y)

  expect_identical(colnames(est(y)[[1]][[1]][[1]]),c('estimate','10%','90%'))
  expect_identical(colnames(est(y,ui_level=.95)[[1]][[1]][[1]]),c('estimate','2.5%','97.5%'))
  expect_identical(y$modelframe,z1$modelframe)

  x <- prepare_data(otu_table=DAT$ABUND,
                    rows_are_taxa=FALSE,
                    tax_table=DAT$TAX,
                    metadata=DAT$META,
                    formula=~s(Day) + Multi,
                    refs='1',
                    cn_normalize=FALSE,
                    drop=TRUE)

  y <- find_topics(x,K=5,init_type='Spectral',tol=1e-03)
  z2 <- est(y)

  expect_true(mean(abs(z1$topic_effects$Multi2$est[,1]-z2$topic_effects$Multi2$est[,1])) < .1)
  expect_true(mean(abs(z1$topic_effects$Multi3$est[,1]-z2$topic_effects$Multi3$est[,1])) < .1)
  expect_true(mean(abs(z1$topic_effects$Day$est[,1]-z2$topic_effects$Day$est[,1])) < .1)

})


context('find topics')

test_that('find_topics returns correct results for different params',{

  skip_on_cran()

  DAT <- readRDS(system.file('testdata','otudata.rds',package='themetagenomics'))

  set.seed(423)

  x <- prepare_data(otu_table=DAT$OTU,
                    rows_are_taxa=FALSE,
                    tax_table=DAT$TAX,
                    metadata=DAT$META,
                    formula=~DIAGNOSIS,
                    refs='Not IBD',
                    cn_normalize=TRUE,
                    drop=TRUE)

  set.seed(1)

  t1 <- find_topics(x,K=5,init_type='Spectral')
  t2 <- find_topics(x,K=5,init_type='LDA')
  t3 <- find_topics(x,K=5,init_type='Random')

  expect_identical(names(t1),names(t2))
  expect_identical(names(t1),names(t3))

  expect_true(cor(sort(colSums(t1$fit$theta)),
                  sort(colSums(t2$fit$theta))) > .5)
  expect_true(cor(sort(colSums(t1$fit$theta)),
                  sort(colSums(t3$fit$theta))) > .5)

  x <- prepare_data(otu_table=DAT$OTU,
                    rows_are_taxa=FALSE,
                    tax_table=DAT$TAX,
                    metadata=DAT$META,
                    cn_normalize=TRUE,
                    drop=TRUE)

  set.seed(1)

  t1 <- find_topics(x,K=5,init_type='Spectral')
  t2 <- find_topics(x,K=5,init_type='LDA')
  t3 <- find_topics(x,K=5,init_type='Random')

  expect_identical(names(t1),names(t2))
  expect_identical(names(t1),names(t3))

  expect_true(cor(sort(colSums(t1$fit$theta)),
                  sort(colSums(t2$fit$theta))) > .5)
  expect_true(cor(sort(colSums(t1$fit$theta)),
                  sort(colSums(t3$fit$theta))) > .5)

})

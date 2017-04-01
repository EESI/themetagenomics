context('find topics')

test_that('find_topics returns correct results for different params',{

  set.seed(423)

  ID1 <- sample(which(GEVERS$META$DIAGNOSIS == 'CD'),10)
  ID2 <- sample(which(GEVERS$META$DIAGNOSIS == 'Not IBD'),10)
  ID <- c(ID1,ID2)
  OTU <- GEVERS$OTU[ID,]
  OTU <- OTU[,colSums(OTU) > 100]
  META <- GEVERS$META[rownames(OTU),]
  TAX <- GEVERS$TAX[colnames(OTU),]

  x <- prepare_data(otu_table=OTU,
                    rows_are_taxa=FALSE,
                    tax_table=TAX,
                    metadata=META,
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
                  sort(colSums(t2$fit$theta))) > .8)
  expect_true(cor(sort(colSums(t1$fit$theta)),
                  sort(colSums(t3$fit$theta))) > .6)

  x <- prepare_data(otu_table=OTU,
                    rows_are_taxa=FALSE,
                    tax_table=TAX,
                    metadata=META,
                    cn_normalize=TRUE,
                    drop=TRUE)

  set.seed(1)

  t1 <- find_topics(x,K=5,init_type='Spectral')
  t2 <- find_topics(x,K=5,init_type='LDA')
  t3 <- find_topics(x,K=5,init_type='Random')

  expect_identical(names(t1),names(t2))
  expect_identical(names(t1),names(t3))

  expect_true(cor(sort(colSums(t1$fit$theta)),
                  sort(colSums(t2$fit$theta))) > .8)
  expect_true(cor(sort(colSums(t1$fit$theta)),
                  sort(colSums(t3$fit$theta))) > .6)

})

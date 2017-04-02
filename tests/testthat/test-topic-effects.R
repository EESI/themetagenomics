context('estimate topic effects')

test_that('est.topics returns correct results for different formulae',{

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

  y <- find_topics(x,K=5,init_type='Spectral')

  z1 <- est(y)

  expect_identical(colnames(est(y)[[1]][[1]][[1]]),c('estimate','10%','90%'))
  expect_identical(colnames(est(y,ui_level=.95)[[1]][[1]][[1]]),c('estimate','2.5%','97.5%'))
  expect_identical(y$modelframe,z1$modelframe)

  expect_error(est(y,metadata=META[1:7,],formula=~PCDAI))

  x <- prepare_data(otu_table=OTU,
                    rows_are_taxa=FALSE,
                    tax_table=TAX,
                    metadata=META,
                    formula=~PCDAI,
                    refs='Not IBD',
                    cn_normalize=TRUE,
                    drop=TRUE)

  y <- find_topics(x,K=5,init_type='Spectral')
  expect_identical(y$modelframe,est(y)$modelframe)

  expect_error(est(y,metadata=META[1:7,],formula=~DIAGNOSIS))
  expect_warning(est(y,metadata=META,formula=~DIAGNOSIS))

  z2 <- est(y,metadata=META,formula=~DIAGNOSIS,refs='Not IBD')

  expect_identical(z1$modelframe,z2$modelframe)
  expect_true(mean(abs(z1$topic_effects$DIAGNOSISCD$est[,1] - z2$topic_effects$DIAGNOSISCD$est[,1])) < .1)

  set.seed(423)

  ID1 <- sample(which(DAVID$META$Site == 'UBERON:feces'),10)
  ID2 <- sample(which(DAVID$META$Site == 'UBERON:saliva'),10)
  ID <- c(ID1,ID2)
  OTU <- DAVID$ABUND[ID,]
  OTU <- OTU[,colSums(OTU) > 100]
  META <- DAVID$META[rownames(OTU),]
  TAX <- DAVID$TAX[colnames(OTU),]

  x <- prepare_data(otu_table=OTU,
                    rows_are_taxa=FALSE,
                    tax_table=TAX,
                    metadata=META,
                    formula=~Site + Day,
                    refs='UBERON:saliva',
                    cn_normalize=FALSE,
                    drop=TRUE)

  y <- find_topics(x,K=5,init_type='Spectral')
  z1 <- est(y)

  expect_identical(colnames(est(y)[[1]][[1]][[1]]),c('estimate','10%','90%'))
  expect_identical(colnames(est(y,ui_level=.95)[[1]][[1]][[1]]),c('estimate','2.5%','97.5%'))
  expect_identical(y$modelframe,z1$modelframe)

  x <- prepare_data(otu_table=OTU,
                    rows_are_taxa=FALSE,
                    tax_table=TAX,
                    metadata=META,
                    formula=~Day + Site,
                    refs='UBERON:saliva',
                    cn_normalize=FALSE,
                    drop=TRUE)

  y <- find_topics(x,K=5,init_type='Spectral')
  z2 <- est(y)

  expect_true(mean(abs(z1$topic_effects$`SiteUBERON:feces`$est[,1]-z2$topic_effects$`SiteUBERON:feces`$est[,1])) < .1)
  expect_true(mean(abs(z1$topic_effects$Day$est[,1]-z2$topic_effects$Day$est[,1])) < .1)

  set.seed(423)

  ID1 <- sample(which(DAVID$META$Site == 'UBERON:feces'),10)
  ID2 <- sample(which(DAVID$META$Site == 'UBERON:saliva'),10)
  ID <- c(ID1,ID2)
  OTU <- DAVID$ABUND[ID,]
  OTU <- OTU[,colSums(OTU) > 100]
  META <- DAVID$META[rownames(OTU),]
  TAX <- DAVID$TAX[colnames(OTU),]

  x <- prepare_data(otu_table=OTU,
                    rows_are_taxa=FALSE,
                    tax_table=TAX,
                    metadata=META,
                    formula=~Site + s(Day),
                    refs='UBERON:saliva',
                    cn_normalize=FALSE,
                    drop=TRUE)

  y <- find_topics(x,K=5,init_type='Spectral')
  z1 <- est(y)

  expect_identical(colnames(est(y)[[1]][[1]][[1]]),c('estimate','10%','90%'))
  expect_identical(colnames(est(y,ui_level=.95)[[1]][[1]][[1]]),c('estimate','2.5%','97.5%'))
  expect_identical(y$modelframe,z1$modelframe)

  x <- prepare_data(otu_table=OTU,
                    rows_are_taxa=FALSE,
                    tax_table=TAX,
                    metadata=META,
                    formula=~s(Day) + Site,
                    refs='UBERON:saliva',
                    cn_normalize=FALSE,
                    drop=TRUE)

  y <- find_topics(x,K=5,init_type='Spectral')
  z2 <- est(y)

  expect_true(mean(abs(z1$topic_effects$`SiteUBERON:feces`$est[,1]-z2$topic_effects$`SiteUBERON:feces`$est[,1])) < .1)
  expect_true(mean(abs(z1$topic_effects$Day$est[,1]-z2$topic_effects$Day$est[,1])) < .1)

  set.seed(123)

  multi <- as.character(as.integer(as.factor(DAVID$META$Site):as.factor(DAVID$META$Donor)))

  ID1 <- sample(which(multi == '1'),10)
  ID2 <- sample(which(multi == '2'),10)
  ID3 <- sample(which(multi == '3'),10)
  ID <- c(ID1,ID2,ID3)
  OTU <- DAVID$ABUND[ID,]
  OTU <- OTU[,colSums(OTU) > 100]
  META <- DAVID$META[rownames(OTU),]
  TAX <- DAVID$TAX[colnames(OTU),]
  META$Multi <- as.character(as.integer(as.factor(META$Site):as.factor(META$Donor)))

  x <- prepare_data(otu_table=OTU,
                    rows_are_taxa=FALSE,
                    tax_table=TAX,
                    metadata=META,
                    formula=~Multi + s(Day),
                    refs='1',
                    cn_normalize=FALSE,
                    drop=TRUE)

  y <- find_topics(x,K=5,init_type='Spectral')
  z1 <- est(y)

  expect_identical(colnames(est(y)[[1]][[1]][[1]]),c('estimate','10%','90%'))
  expect_identical(colnames(est(y,ui_level=.95)[[1]][[1]][[1]]),c('estimate','2.5%','97.5%'))
  expect_identical(y$modelframe,z1$modelframe)

  x <- prepare_data(otu_table=OTU,
                    rows_are_taxa=FALSE,
                    tax_table=TAX,
                    metadata=META,
                    formula=~s(Day) + Multi,
                    refs='1',
                    cn_normalize=FALSE,
                    drop=TRUE)

  y <- find_topics(x,K=5,init_type='Spectral')
  z2 <- est(y)

  expect_true(mean(abs(z1$topic_effects$Multi2$est[,1]-z2$topic_effects$Multi2$est[,1])) < .1)
  expect_true(mean(abs(z1$topic_effects$Multi3$est[,1]-z2$topic_effects$Multi3$est[,1])) < .1)
  expect_true(mean(abs(z1$topic_effects$Day$est[,1]-z2$topic_effects$Day$est[,1])) < .1)

})


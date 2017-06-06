context('prediction functions')

test_that('picrust for cog and kegg both work',{

  skip_on_cran()

  DAT <- readRDS(system.file('testdata','otudata.rds',package='themetagenomics'))

  x <- prepare_data(otu_table=DAT$OTU,
                    rows_are_taxa=FALSE,
                    tax_table=DAT$TAX,
                    metadata=DAT$META,
                    formula=~DIAGNOSIS,
                    refs='Not IBD',
                    cn_normalize=TRUE,
                    drop=TRUE)

  y <- find_topics(x,K=5,verbose=FALSE)

  tmp <- tempdir()
  download_ref(tmp,reference='gg_ko',overwrite=FALSE)
  expect_warning(z1 <- predict(y,reference='gg_ko',reference_path=tmp,cn_normalize=TRUE))

  download_ref(tmp,reference='gg_cog',overwrite=FALSE)
  z2 <- predict(y,reference='gg_cog',reference_path=tmp)

  expect_identical(names(z1),names(z2))
  expect_identical(nrow(z1$fxn_table),nrow(z2$fxn_table))
  expect_identical(min(z1$fxn_table),min(z2$fxn_table))
  expect_identical(max(z1$fxn_table) > 1000,max(z2$fxn_table) > 1000)
  expect_identical(min(z1$method_meta) >= 0,min(z2$method_meta) >= 0)
  expect_identical(max(z1$method_meta) <= 1,max(z2$method_meta) <= 1)

})

test_that('t4f works for different params',{

  DAT <- readRDS(system.file('testdata','seqdata.rds',package='themetagenomics'))

  x <- prepare_data(otu_table=DAT$ABUND,
                    rows_are_taxa=FALSE,
                    tax_table=DAT$TAX,
                    metadata=DAT$META,
                    formula=~Day,
                    cn_normalize=FALSE,
                    drop=TRUE)

  y <- find_topics(x,K=5)

  tmp <- tempdir()
  download_ref(tmp,reference='silva_ko',overwrite=FALSE)
  expect_warning(z1 <- predict(y,reference='silva_ko',reference_path=tmp,type='uproc',short=TRUE,cn_normalize=FALSE,sample_normalize=FALSE))
  z2 <- predict(y,reference='silva_ko',reference_path=tmp,type='uproc',short=TRUE,cn_normalize=TRUE,sample_normalize=TRUE)
  z3 <- predict(y,reference='silva_ko',reference_path=tmp,type='uproc',short=FALSE,cn_normalize=TRUE,sample_normalize=FALSE)
  z4 <- predict(y,reference='silva_ko',reference_path=tmp,type='pauda',short=TRUE,cn_normalize=TRUE,sample_normalize=FALSE)

  expect_true(min(z1$method_meta) >= 0 & max(z1$method_meta) <= 1)
  expect_true(min(z2$method_meta) >= 0 & max(z2$method_meta) <= 1)
  expect_true(min(z3$method_meta) >= 0 & max(z3$method_meta) <= 1)
  expect_true(min(z4$method_meta) >= 0 & max(z4$method_meta) <= 1)

  expect_true(all(c(max(z1$fxn_table),max(z2$fxn_table),max(z3$fxn_table),max(z4$fxn_table)) == 100))
  expect_identical(names(z1),names(z2))
  expect_identical(names(z2),names(z3))
  expect_identical(names(z3),names(z4))

  overlap <- intersect(intersect(intersect(colnames(z1$fxn_table),colnames(z2$fxn_table)),colnames(z3$fxn_table)),colnames(z4$fxn_table))
  expect_true(length(overlap) > 1000)

  expect_true(all(length(z1$fxn_meta[overlap]) == c(length(z2$fxn_meta[overlap]),length(z3$fxn_meta[overlap]),length(z4$fxn_meta[overlap]))))

})

context("prepare data")

test_that('prepare_data works with no formula',{

  DAT <- readRDS(system.file('testdata','otudata.rds',package='themetagenomics'))

  x <- prepare_data(otu_table=DAT$OTU,
                    rows_are_taxa=FALSE,
                    cn_normalize=TRUE,
                    drop=TRUE)
  expect_identical(list('matrix',NULL,NULL,NULL,NULL,NULL,NULL),c(class(x[[1]])[1],unname(x[2:7])))

  x <- prepare_data(otu_table=DAVID$ABUND,
                    rows_are_taxa=FALSE,
                    cn_normalize=FALSE,
                    drop=TRUE)
  expect_identical(list('matrix',NULL,NULL,NULL,NULL,NULL,NULL),c(class(x[[1]])[1],unname(x[2:7])))

  expect_warning(prepare_data(otu_table=DAVID$ABUND,
                              rows_are_taxa=FALSE,
                              cn_normalize=TRUE,
                              drop=TRUE))
})

test_that('prepare_data works with binary factor',{

  DAT <- readRDS(system.file('testdata','otudata.rds',package='themetagenomics'))

  x <- prepare_data(otu_table=DAT$OTU,
                    rows_are_taxa=FALSE,
                    tax_table=DAT$TAX,
                    metadata=DAT$META,
                    formula=~DIAGNOSIS,
                    refs='Not IBD',
                    cn_normalize=TRUE,
                    drop=TRUE)
  # expect_identical(list(c('matrix','array'),c('matrix','array'),
  #                       'data.frame','formula','character','NULL','data.frame','list'),
  #                  unname(lapply(x,class)))
  expect_true(colnames(x$metadata) != colnames(x$modelframe))

})

test_that('prepare_data works with continuous covariate',{

  DAT <- readRDS(system.file('testdata','otudata.rds',package='themetagenomics'))
  DAT$META$PCDAI[1] <- NA

  x <- prepare_data(otu_table=DAT$OTU,
                    rows_are_taxa=FALSE,
                    tax_table=DAT$TAX,
                    metadata=DAT$META,
                    formula=~PCDAI,
                    cn_normalize=TRUE,
                    drop=TRUE)
  # expect_identical(list(c('matrix','array'),c('matrix','array'),
  #                       'data.frame','formula','NULL','data.frame','list'),
  #                  unname(lapply(x,class)))
  expect_true(any(is.na(DAT$META$PCDAI)))
  expect_false(any(is.na(x$metadata$PCDAI)))

})

test_that('prepare_data works with spline covariate',{

  DAT <- readRDS(system.file('testdata','otudata.rds',package='themetagenomics'))

  x <- prepare_data(otu_table=DAT$OTU,
                    rows_are_taxa=FALSE,
                    tax_table=DAT$TAX,
                    metadata=DAT$META,
                    formula=~s(PCDAI),
                    cn_normalize=TRUE,
                    drop=FALSE)
  # expect_identical(list(c('matrix','array'),c('matrix','array'),
  #                       'data.frame','formula','list','data.frame','list'),
  #                  unname(lapply(x,class)))
  expect_equal(list('PCDAI',
                    's',
                    model.frame(~s(PCDAI),na.omit(DAT$META)),
                    ~PCDAI),
               list(x$splineinfo[[1]][[1]][[1]],
                    x$splineinfo[[1]][[1]][[2]],
                    x$splineinfo[[1]][[1]][[3]],
                    x$splineinfo[[2]]))
})

test_that('prepare_data works with multiclass factor',{

  DAT <- readRDS(system.file('testdata','seqdata.rds',package='themetagenomics'))

  # expect_warning(prepare_data(otu_table=DAT$ABUND,
  #                             rows_are_taxa=FALSE,
  #                             tax_table=DAT$TAX,
  #                             metadata=DAT$META,
  #                             formula=~Multi,
  #                             cn_normalize=FALSE,
  #                             drop=TRUE))

  x <- prepare_data(otu_table=DAT$ABUND,
                    rows_are_taxa=FALSE,
                    tax_table=DAT$TAX,
                    metadata=DAT$META,
                    refs='2',
                    formula=~Multi,
                    cn_normalize=FALSE,
                    drop=TRUE)
  # expect_identical(list(c('matrix','array'),c('matrix','array'),
  #                       'data.frame','formula','character','NULL','data.frame','list'),
  #                  unname(lapply(x,class)))
  expect_true(ncol(x$modelframe) > ncol(x$metadata))

})

test_that('prepare_data works with multiclass factor and continuous',{

  DAT <- readRDS(system.file('testdata','seqdata.rds',package='themetagenomics'))

  x <- prepare_data(otu_table=DAT$ABUND,
                    rows_are_taxa=FALSE,
                    tax_table=DAT$TAX,
                    metadata=DAT$META,
                    refs='2',
                    formula=~Multi + Day,
                    cn_normalize=FALSE,
                    drop=TRUE)
  # expect_identical(list(c('matrix','array'),c('matrix','array'),
  #                       'data.frame','formula','character','NULL','data.frame','list'),
  #                  unname(lapply(x,class)))
  expect_true(ncol(x$modelframe) > ncol(x$metadata))
  y <- prepare_data(otu_table=DAT$ABUND,
                    rows_are_taxa=FALSE,
                    tax_table=DAT$TAX,
                    metadata=DAT$META,
                    refs='2',
                    formula=~Day + Multi,
                    cn_normalize=FALSE,
                    drop=TRUE)
  expect_equal(x$modelframe[,sort(colnames(x$modelframe))],
               y$modelframe[,sort(colnames(y$modelframe))])

})

test_that('prepare_data works with multiclass factor and spline',{

  DAT <- readRDS(system.file('testdata','seqdata.rds',package='themetagenomics'))

  x <- prepare_data(otu_table=DAT$ABUND,
                    rows_are_taxa=FALSE,
                    tax_table=DAT$TAX,
                    metadata=DAT$META,
                    refs='2',
                    formula=~Multi + s(Day),
                    cn_normalize=FALSE,
                    drop=TRUE)
  # expect_identical(list(c('matrix','array'),c('matrix','array'),
  #                       'data.frame','formula','character','list','data.frame','list'),
  #                  unname(lapply(x,class)))
  expect_true(ncol(x$modelframe) > ncol(x$metadata))
  y <- prepare_data(otu_table=DAT$ABUND,
                    rows_are_taxa=FALSE,
                    tax_table=DAT$TAX,
                    metadata=DAT$META,
                    refs='2',
                    formula=~s(Day) + Multi,
                    cn_normalize=FALSE,
                    drop=TRUE)
  expect_equal(x$modelframe[,sort(colnames(x$modelframe))],
               y$modelframe[,sort(colnames(y$modelframe))])
  expect_equal(list('Day',
                    's',
                    model.frame(~s(Day),na.omit(DAT$META)),
                    ~Multi + Day),
               list(x$splineinfo[[1]][[1]][[1]],
                    x$splineinfo[[1]][[1]][[2]],
                    x$splineinfo[[1]][[1]][[3]],
                    x$splineinfo[[2]]))
  expect_equal(list('Day',
                    's',
                    model.frame(~s(Day),na.omit(DAT$META)),
                    ~Day + Multi),
               list(y$splineinfo[[1]][[1]][[1]],
                    y$splineinfo[[1]][[1]][[2]],
                    y$splineinfo[[1]][[1]][[3]],
                    y$splineinfo[[2]]))
})

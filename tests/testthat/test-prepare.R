context("prepare data")

test_that('prepare_data works with no formula',{
  x <- prepare_data(otu_table=GEVERS$OTU,
                    rows_are_taxa=FALSE,
                    cn_normalize=TRUE,
                    drop=TRUE)
  expect_identical(list('matrix',NULL,NULL,NULL,NULL,NULL,NULL),c(class(x[[1]]),unname(x[2:7])))

  x <- prepare_data(otu_table=DAVID$ABUND,
                    rows_are_taxa=FALSE,
                    cn_normalize=FALSE,
                    drop=TRUE)
  expect_identical(list('matrix',NULL,NULL,NULL,NULL,NULL,NULL),c(class(x[[1]]),unname(x[2:7])))

  expect_warning(prepare_data(otu_table=DAVID$ABUND,
                              rows_are_taxa=FALSE,
                              cn_normalize=TRUE,
                              drop=TRUE))
})

test_that('prepare_data works with binary factor',{
  x <- prepare_data(otu_table=GEVERS$OTU,
                    rows_are_taxa=FALSE,
                    tax_table=GEVERS$TAX,
                    metadata=GEVERS$META,
                    formula=~DIAGNOSIS,
                    refs='Not IBD',
                    cn_normalize=TRUE,
                    drop=TRUE)
  expect_identical(list('matrix','matrix','data.frame','formula','character','NULL','data.frame'),
                   unname(lapply(x,class)))
  expect_true(colnames(x$metadata) != colnames(x$modelframe))
})

test_that('prepare_data works with continuous covariate',{
  x <- prepare_data(otu_table=GEVERS$OTU,
                    rows_are_taxa=FALSE,
                    tax_table=GEVERS$TAX,
                    metadata=GEVERS$META,
                    formula=~PCDAI,
                    cn_normalize=TRUE,
                    drop=TRUE)
  expect_identical(list('matrix','matrix','data.frame','formula','NULL','data.frame'),
                   unname(lapply(x,class)))
  expect_true(any(is.na(GEVERS$META$PCDAI)))
  expect_false(any(is.na(x$metadata$PCDAI)))
})

test_that('prepare_data works with spline covariate',{
  x <- prepare_data(otu_table=GEVERS$OTU,
                    rows_are_taxa=FALSE,
                    tax_table=GEVERS$TAX,
                    metadata=GEVERS$META,
                    formula=~s(PCDAI),
                    cn_normalize=TRUE,
                    drop=TRUE)
  expect_identical(list('matrix','matrix','data.frame','formula','list','data.frame'),
                   unname(lapply(x,class)))
  expect_equal(list('PCDAI',
                    's',
                    model.frame(~s(PCDAI),na.omit(GEVERS$META)),
                    ~PCDAI),
               list(x$splineinfo[[1]][[1]][[1]],
                     x$splineinfo[[1]][[1]][[2]],
                     x$splineinfo[[1]][[1]][[3]],
                     x$splineinfo[[2]]))
})

test_that('prepare_data works with multiclass factor',{
  DAVID$META$SiteDonor <- with(DAVID$META,as.factor(Site):as.factor(Donor))
  expect_warning(prepare_data(otu_table=DAVID$ABUND,
                              rows_are_taxa=FALSE,
                              tax_table=DAVID$TAX,
                              metadata=DAVID$META,
                              formula=~SiteDonor,
                              cn_normalize=FALSE,
                              drop=TRUE))
  DAVID$META$SiteDonor <- as.character(as.integer(DAVID$META$SiteDonor))
  x <- prepare_data(otu_table=DAVID$ABUND,
                    rows_are_taxa=FALSE,
                    tax_table=DAVID$TAX,
                    metadata=DAVID$META,
                    refs='2',
                    formula=~SiteDonor,
                    cn_normalize=FALSE,
                    drop=TRUE)
  expect_identical(list('matrix','matrix','data.frame','formula','character','NULL','data.frame'),
                   unname(lapply(x,class)))
  expect_true(ncol(x$modelframe) > ncol(x$metadata))
})

test_that('prepare_data works with multiclass factor and continuous',{
  DAVID$META$SiteDonor <- with(DAVID$META,as.factor(Site):as.factor(Donor))
  DAVID$META$SiteDonor <- as.character(as.integer(DAVID$META$SiteDonor))
  x <- prepare_data(otu_table=DAVID$ABUND,
                    rows_are_taxa=FALSE,
                    tax_table=DAVID$TAX,
                    metadata=DAVID$META,
                    refs='2',
                    formula=~SiteDonor + Day,
                    cn_normalize=FALSE,
                    drop=TRUE)
  expect_identical(list('matrix','matrix','data.frame','formula','character','NULL','data.frame'),
                   unname(lapply(x,class)))
  expect_true(ncol(x$modelframe) > ncol(x$metadata))
  y <- prepare_data(otu_table=DAVID$ABUND,
                    rows_are_taxa=FALSE,
                    tax_table=DAVID$TAX,
                    metadata=DAVID$META,
                    refs='2',
                    formula=~Day + SiteDonor,
                    cn_normalize=FALSE,
                    drop=TRUE)
  expect_equal(x$modelframe[,c('SiteDonor1','SiteDonor3','Day')],
               y$modelframe[,c('SiteDonor1','SiteDonor3','Day')])
})

test_that('prepare_data works with multiclass factor and spline',{
  DAVID$META$SiteDonor <- with(DAVID$META,as.factor(Site):as.factor(Donor))
  DAVID$META$SiteDonor <- as.character(as.integer(DAVID$META$SiteDonor))
  x <- prepare_data(otu_table=DAVID$ABUND,
                    rows_are_taxa=FALSE,
                    tax_table=DAVID$TAX,
                    metadata=DAVID$META,
                    refs='2',
                    formula=~SiteDonor + s(Day),
                    cn_normalize=FALSE,
                    drop=TRUE)
  expect_identical(list('matrix','matrix','data.frame','formula','character','list','data.frame'),
                   unname(lapply(x,class)))
  expect_true(ncol(x$modelframe) > ncol(x$metadata))
  y <- prepare_data(otu_table=DAVID$ABUND,
                    rows_are_taxa=FALSE,
                    tax_table=DAVID$TAX,
                    metadata=DAVID$META,
                    refs='2',
                    formula=~s(Day) + SiteDonor,
                    cn_normalize=FALSE,
                    drop=TRUE)
  expect_equal(x$modelframe[,c('SiteDonor1','SiteDonor3','Day')],
               y$modelframe[,c('SiteDonor1','SiteDonor3','Day')])
  expect_equal(list('Day',
                    's',
                    model.frame(~s(Day),na.omit(DAVID$META)),
                    ~SiteDonor + Day),
               list(x$splineinfo[[1]][[1]][[1]],
                    x$splineinfo[[1]][[1]][[2]],
                    x$splineinfo[[1]][[1]][[3]],
                    x$splineinfo[[2]]))
  expect_equal(list('Day',
                    's',
                    model.frame(~s(Day),na.omit(DAVID$META)),
                    ~Day + SiteDonor),
               list(y$splineinfo[[1]][[1]][[1]],
                    y$splineinfo[[1]][[1]][[2]],
                    y$splineinfo[[1]][[1]][[3]],
                    y$splineinfo[[2]]))
})

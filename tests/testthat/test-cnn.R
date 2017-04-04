context("copy number normalize")

test_that('cnn returns matrix of same dimensions',{

  DAT <- readRDS(system.file('testdata','otudata.rds',package='themetagenomics'))

  expect_equal(cnn(DAT$OTU,rows_are_taxa=FALSE,drop=FALSE),
               cnn(t(DAT$OTU),rows_are_taxa=TRUE,drop=FALSE))

    expect_equal(cnn(DAT$OTU,rows_are_taxa=FALSE,drop=TRUE),
               cnn(t(DAT$OTU),rows_are_taxa=TRUE,drop=TRUE))
})

test_that('cnn returns warning and same matrix with naming issue',{

  DAT <- readRDS(system.file('testdata','otudata.rds',package='themetagenomics'))

  expect_warning(cnn(unname(DAT$OTU),rows_are_taxa=FALSE,drop=FALSE))
  expect_equal(unname(DAT$OTU),
               suppressWarnings(cnn(unname(DAT$OTU),
                                    rows_are_taxa=FALSE,drop=FALSE)))

  DAT <- readRDS(system.file('testdata','seqdata.rds',package='themetagenomics'))

  expect_warning(cnn(unname(DAT$ABUND),rows_are_taxa=FALSE,drop=FALSE))
  expect_equal(unname(DAT$ABUND),
               suppressWarnings(cnn(unname(DAT$ABUND),
                                    rows_are_taxa=FALSE,drop=FALSE)))
})

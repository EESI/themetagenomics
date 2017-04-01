context("copy number normalize")

test_that('cnn returns matrix of same dimensions',{
  expect_equal(cnn(GEVERS$OTU,rows_are_taxa=FALSE,drop=FALSE),
               cnn(t(GEVERS$OTU),rows_are_taxa=TRUE,drop=FALSE))

    expect_equal(cnn(GEVERS$OTU,rows_are_taxa=FALSE,drop=TRUE),
               cnn(t(GEVERS$OTU),rows_are_taxa=TRUE,drop=TRUE))
})

test_that('cnn returns warning and same matrix with naming issue',{
  expect_warning(cnn(unname(GEVERS$OTU),rows_are_taxa=FALSE,drop=FALSE))
  expect_equal(unname(GEVERS$OTU),
               suppressWarnings(cnn(unname(GEVERS$OTU),
                                    rows_are_taxa=FALSE,drop=FALSE)))
  expect_warning(cnn(unname(DAVID$ABUND),rows_are_taxa=FALSE,drop=FALSE))
  expect_equal(unname(DAVID$ABUND),
               suppressWarnings(cnn(unname(DAVID$ABUND),
                                    rows_are_taxa=FALSE,drop=FALSE)))
})

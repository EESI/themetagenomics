context('misc functions')

test_that('pretty_taxa_names works for all levels',{

  DAT <- readRDS(system.file('testdata','otudata.rds',package='themetagenomics'))

  expect_identical(length(pretty_taxa_names(DAT$TAX)),nrow(DAT$TAX))
  expect_false(any(is.na(pretty_taxa_names(DAT$TAX))))
  expect_equal(unique(sapply(strsplit(pretty_taxa_names(DAT$TAX),' '),length)),2)

})

test_that('renaming taxa to other always returns other',{

  DAT <- readRDS(system.file('testdata','otudata.rds',package='themetagenomics'))
  OTHER <- rename_taxa_to_other(DAT$OTU,DAT$TAX)

  expect_equal(unique(sapply(gregexpr('\\W+', unlist(OTHER[,-c(1,7)])),length)),1)

})

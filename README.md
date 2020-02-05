
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Build
Status](https://travis-ci.com/sw1/themetagenomics.svg?token=8r1TnJBy8TyidNrmbPpa&branch=master)](https://travis-ci.com/sw1/themetagenomics)
[![codecov](https://codecov.io/gh/sw1/themetagenomics/branch/master/graph/badge.svg?token=pmjXMfuHrw)](https://codecov.io/gh/sw1/themetagenomics)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/themetagenomics)](https://cran.r-project.org/package=themetagenomics)

## Overview

themetagenomics provides functions to explore topics generated from 16S
rRNA sequencing information on both the abundance and functional levels.
It also provides an R implementation of PICRUSt and wraps Tax4fun,
giving users a choice for their functional prediction strategy.

## Relevant Software and Literature

[Stephen Woloszynek, Joshua Chang Mell, Zhengqiao Zhao, Gideon Simpson,
Michael P. O’Connor, and Gail Rosen. Exploring thematic structure and
predicted functionality of 16S rRNA amplicon data. 2019. PLoS
ONE 14(12); doi:
https://doi.org/10.1371/journal.pone.0219235](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0219235)

[Margaret E. Roberts, Brandon M. Stewart, and Dustin Tingley (2017).
stm: R Package for Structural Topic
Models, 2016.](http://www.structuraltopicmodel.com)

[Morgan Langille et al. Predictive functional profiling of microbial
communities using 16S rRNA marker gene sequences.
PICRUSt 1.1.1](http://picrust.github.io/picrust/)

[Kathrin P. Aßhauer, Bernd Wemheuer, Rolf Daniel, and Peter Meinicke.
Tax4Fun: predicting functional profiles from metagenomic 16S rRNA
data](http://tax4fun.gobics.de/)

[Stan Development Team. 2016. RStan: the R interface to Stan. R package
version 2.14.1., 2017.](http://mc-stan.org)

## Installation

``` r
# To install the release version of Themetagenomics:
install.packages('themetagenomics')

# To install the developmental version via Github:
# install.packages('devtools')
devtools::install_github('EESI/themetagenomics',build_vignettes=TRUE)
```

## Development

For future feature requests or suggestions, request an issue:

<https://github.com/EESI/themetagenomics/issues>

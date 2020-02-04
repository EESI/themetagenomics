
<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Build Status](https://travis-ci.com/sw1/themetagenomics.svg?token=8r1TnJBy8TyidNrmbPpa&branch=master)](https://travis-ci.com/sw1/themetagenomics) [![codecov](https://codecov.io/gh/sw1/themetagenomics/branch/master/graph/badge.svg)](https://codecov.io/gh/sw1/themetagenomics) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/themetagenomics)](https://cran.r-project.org/package=themetagenomics)

Overview
--------

themetagenomics provides functions to explore topics generated from 16S rRNA sequencing information on both the abundance and functional levels. It also provides an R implementation of PICRUSt and wraps Tax4fun, giving users a choice for their functional prediction strategy.

Relevant Software and Literature
--------------------------------

[Stephen Woloszynek, Joshua Chang Mell, Gideon Simpson, and Gail Rosen. Exploring thematic structure in 16S rRNA marker gene surveys. 2017. bioRxiv 146126; doi: https://doi.org/10.1101/146126](http://biorxiv.org/content/early/2017/06/05/146126.article-info)

Stephen Woloszynek, Zhengqiao Zhao, Gideon Simpson, Joshua Chang Mell, and Gail Rosen. Gauging the use of topic models as a means to understand microbiome data structure. TBD, 2017

[Margaret E. Roberts, Brandon M. Stewart, and Dustin Tingley (2017). stm: R Package for Structural Topic Models, 2016.](http://www.structuraltopicmodel.com)

[Morgan Langille et al. Predictive functional profiling of microbial communities using 16S rRNA marker gene sequences. PICRUSt 1.1.1](http://picrust.github.io/picrust/)

[Kathrin P. Aßhauer, Bernd Wemheuer, Rolf Daniel, and Peter Meinicke. Tax4Fun: predicting functional profiles from metagenomic 16S rRNA data](http://tax4fun.gobics.de/)

[Stan Development Team. 2016. RStan: the R interface to Stan. R package version 2.14.1., 2017.](http://mc-stan.org)

Installation
------------

``` r
# To install the release version of Themetagenomics:
install.packages('themetagenomics')

# To install the developmental version via Github:
# install.packages('devtools')
devtools::install_github('EESI/themetagenomics',build_vignettes=TRUE)
```

Development
-----------

For future feature requests or suggestions, request an issue:

<https://github.com/EESI/themetagenomics/issues>

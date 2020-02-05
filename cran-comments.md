This is a submission (first release)

---

## Test environments
* OSX (travis-ci): R-release (3.6.2) [2/4/2020]
* ubuntu (travis-ci): R-release (3.6.2) [2/4/2020]
* win (travis-ci): R-release (3.6.2) [2/4/2020]
* ubuntu 16.10 (local machine): R 3.6.2
* win-builder: R-release (3.6.1), R-devel (4.0.0) [2/4/2020]

## R CMD check results

0 errors | 0 warnings | 0 notes

## License 

* License components with restrictions and base license permitting such:
  MIT + file LICENSE
  
## Downstream dependencies

* None (first release)

## Resubmission 1

* URLs: removed empty URLS; change biorxiv doi url to direct url; changed citation doi to biorxiv url
* Native routines and symbol search: ran package_native_routine_registration_skeleton() in R 3.4 and
  placed output in src/themetagenomics_init.c; added .registration = TRUE to useDynLib() in NAMESPACE
* Examples: wrapped all examples exceeding 10s with \dontrun
* Tests: added skip_on_cran() to all long tests

## Resubmission 2

* Description: removed all mentions of the word 'library'
* Title: changed 16S to 16s per toTitleCase(); note that the term 16S should be capitalized, see wikipedia
  for the term 16S rRNA, but I changed it to 16s for this resubmission.

## Resubmission 3

* Adjusted any use of class(x) to prevent issues from length 2 (matrix, array) output if matrix.


<!-- README.md is generated from README.Rmd. Please edit that file -->

# metaboprep2

<!-- badges: start -->

[![Codecov test
coverage](https://codecov.io/gh/nicksunderland/metaboprep2/graph/badge.svg)](https://app.codecov.io/gh/nicksunderland/metaboprep2)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

`metaboprep2` is an R package designed for processing and managing
metabolomic data efficiently. It provides structured storage and access
to metabolomic datasets, including sample metadata, feature metadata,
and metabolite intensity data. Whilst the class structure is new, the
underlying routines are based on the excellent work from the authors of
[metaboprep](https://github.com/MRCIEU/metaboprep). Until this project
is stable and tested it is recommended you use the original `metaboprep`
pipeline.

## Installation

You can install the development version of metaboprep2 from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("nicksunderland/metaboprep2")

# or 
remotes::install_github("nicksunderland/metaboprep2")
```

## Getting Started

Please see the dedicated
[vignettes](https://nicksunderland.github.io/metaboprep2/) for common
use scenarios.

## TODO

- (CLI)\[<https://github.com/r-lib/cli>\] - nice command line output for
  processing.  
- HTML report output
- bug where you cant rerun QC on the same layer due to sparse matrix
  error
- question - why does the report have 3 ‘featuremis\_%NUM’ sections that
  do the same thing?
- question - same for chunks like ‘sample_missingness_3’

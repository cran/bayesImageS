# bayesImageS 0.4-1

## Changes

* removed the deprecated `slices` parameter from `mcmcPotts()`

* some minor edits to the `Roxygen` documentation

* added ORCID to DESCRIPTION

* removed the R package `PottsUtils` from Suggested packages; added `mcmcse`

* moved `README-chunkname.png` images from root directory to `inst/images`

# bayesImageS 0.4-0

## New features

* Functions getNeighbors(), getEdges(), getBlocks() contributed by Dai Feng. These are currently exact copies of the equivalent functions in the R package PottsUtils. The reason for duplicating this code is to avoid picking up dependencies on miscF, R2jags, BRugs, etc.
    + <https://cran.r-project.org/package=PottsUtils>

* Added new option to swNoData and mcmcPottsNoData to initialise the labels using either random or deterministic values. This can be useful when comparing the convergence of the Swendsen-Wang and Gibbs sampling algorithms to known distributions (e.g. when $\beta$=0 or $\infty$)

* Additional documentation in README.md with an example of usage

## Bug fixes

* Fixed compile errors on SPARC Solaris
    + <https://www.r-project.org/nosvn/R.check/r-patched-solaris-x86/bayesImageS-00check.html>

* Fixed compile warnings from icpc (Intel Parallel Studio XE) on Linux

* Added bayesImageS_init.cpp to call R_registerRoutines and R_useDynamicSymbols

# bayesImageS 0.3-3

* First version released on CRAN

# bayesImageS 0.1-21

* First public release on RunMyCode.org:
    + <http://www.runmycode.org/companion/view/546>

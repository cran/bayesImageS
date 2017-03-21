
# bayesImageS 0.4-0

## New features

* Functions getNeighbors(), getEdges(), getBlocks() contributed by Dai Feng. These are currently exact copies of the equivalent functions in the R package PottsUtils. The reason for duplicating this code is to avoid picking up dependencies on miscF, R2jags, BRugs, etc.
    + https://cran.r-project.org/package=PottsUtils

* Added new option to swNoData and mcmcPottsNoData to initialise the labels using either random or deterministic values. This can be useful when comparing the convergence of the Swendsen-Wang and Gibbs sampling algorithms to known distributions (e.g. when beta=0 or infinity)

* Additional documentation in README.md with an example of usage

## Bug fixes

* attempted to fix compile errors on SPARC Solaris
    + https://www.r-project.org/nosvn/R.check/r-patched-solaris-sparc/bayesImageS-00install.html

* fixed compile warnings from icpc (Intel Parallel Studio XE)

* added bayesImageS_init.cpp to call R_registerRoutines and R_useDynamicSymbols

# bayesImageS 0.3-3

* first version released on CRAN

# bayesImageS 0.1-21

* first public release on RunMyCode.org:
    + http://www.runmycode.org/companion/view/546
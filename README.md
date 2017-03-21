
<!--- README.md is generated from README.Rmd. Please edit that file -->
Various algorithms for segmentation of 2D and 3D images, such as computed tomography (CT) and satellite remote sensing. This R package implements Bayesian image analysis using the hidden Potts/Ising model with external field prior. Latent labels are sampled using chequerboard Gibbs sampling or Swendsen-Wang. Algorithms for the smoothing parameter include pseudolikelihood, path sampling, the exchange algorithm, and approximate Bayesian computation (ABC-MCMC and ABC-SMC).

For example, to generate synthetic data for a known value of *β*:

``` r
set.seed(123456)
library(bayesImageS)
neigh <- matrix(c(5,2,5,3,  1,5,5,4,  5,4,1,5,  3,5,2,5), nrow=4, ncol=4, byrow=TRUE)
blocks <- list(c(1,4), c(2,3))
res.sw <- swNoData(0.5, 3, neigh, blocks, niter=200)
```

Now add some Gaussian noise to the labels, according to the prior:

``` r
priors <- list()
priors$k <- 3
priors$mu <- c(-2, 0, 2)
priors$mu.sd <- rep(0.5,priors$k)
priors$sigma <- rep(0.25,priors$k)
priors$sigma.nu <- rep(3, priors$k)
priors$beta <- c(0,1.3)

m0 <- sort(rnorm(priors$k,priors$mu,priors$mu.sd))
SS0 <- priors$sigma.nu*priors$sigma^2
s0 <- 1/sqrt(rgamma(priors$k,priors$sigma.nu/2,SS0/2))
z <- max.col(res.sw$z)[1:nrow(neigh)]
y <- m0[z] + rnorm(nrow(neigh),0,s0[z])
```

Image segmentation using SMC-ABC:

``` r
res.smc <- smcPotts(y, neigh, blocks, priors=priors)
#> Initialization took 1sec
#> Iteration 1
#> previous epsilon 2 and ESS 10000 (target: 9500)
#> Took 8sec to update epsilon=9.33264e-302 (ESS=9796.24)
#> Took 1sec for 9105 RWMH updates (bw=0.519431)
#> Took 1sec for 10000 iterations to calculate S(z)=2
# pixel classifications
pred <- res.smc$alloc/rowSums(res.smc$alloc)
seg <- max.col(res.smc$alloc) # posterior mode (0-1 loss)
all.equal(seg, z)
#> [1] TRUE
mean(res.smc$beta)
#> [1] 0.649113
apply(res.smc$mu, 2, range)
#>            [,1]      [,2]       [,3]
#> [1,] -3.2701475 -1.254859 0.08924757
#> [2,] -0.5336341  1.388502 4.04662164
m0
#> [1] -1.14599925 -0.09531646  1.81442878
```

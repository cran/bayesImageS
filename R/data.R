#' Simulation from the Potts model using Swendsen-Wang.
#' 
#' Simulations for a \eqn{100 \times 100} lattice for fixed values
#' of the inverse temperature parameter, \eqn{\beta}.
#' 
#' @format A \code{list} containing 3 variables:
#' \describe{
#'   \item{beta}{values of beta}
#'   \item{matu}{corresponding simulations of the sufficient statistic}
#'   \item{tm}{time taken for these simulations}
#' }
#' @seealso \code{\link{swNoData}}
"simSW"

#' Simulation from the Potts model using Swendsen-Wang.
#' 
#' Simulations for a \eqn{500 \times 500} lattice for fixed values
#' of the inverse temperature parameter, \eqn{\beta}.
#' 
#' @format A \code{list} containing 5 variables:
#' \describe{
#'   \item{0.22}{simulations for \eqn{\beta = 0.22}}
#'   \item{0.44}{simulations for \eqn{\beta = 0.44}}
#'   \item{0.88}{simulations for \eqn{\beta = 0.88}}
#'   \item{1.32}{simulations for \eqn{\beta = 1.32}}
#'   \item{tm}{time taken by the simulations}
#' }
#' @seealso \code{\link{swNoData}}
"synth"

#' Simulation from the Potts model using single-site Gibbs updates.
#' 
#' 100 iterations of Gibbs sampling for a \eqn{500 \times 500} lattice
#' with \eqn{\beta=0.22} and \eqn{k=2}.
#' @format A \code{list} containing 7 variables.
#' @seealso \code{\link{mcmcPotts}}
"res"

#' Simulation from the Potts model using single-site Gibbs updates.
#' 
#' 100 iterations of Gibbs sampling for a \eqn{500 \times 500} lattice
#' with \eqn{\beta=0.44} and \eqn{k=2}.
#' @format A \code{list} containing 7 variables.
#' @seealso \code{\link{mcmcPotts}}
"res2"

#' Simulation from the Potts model using single-site Gibbs updates.
#' 
#' 100 iterations of Gibbs sampling for a \eqn{500 \times 500} lattice
#' with \eqn{\beta=0.88} and \eqn{k=2}.
#' @format A \code{list} containing 7 variables.
#' @seealso \code{\link{mcmcPotts}}
"res3"

#' Simulation from the Potts model using single-site Gibbs updates.
#' 
#' 100 iterations of Gibbs sampling for a \eqn{500 \times 500} lattice
#' with \eqn{\beta=1.32} and \eqn{k=2}.
#' @format A \code{list} containing 7 variables.
#' @seealso \code{\link{mcmcPotts}}
"res4"

#' Simulation from the Potts model using single-site Gibbs updates.
#' 
#' 5000 iterations of Gibbs sampling for a \eqn{500 \times 500} lattice
#' with \eqn{\beta=1.32} and \eqn{k=2}.
#' @format A \code{list} containing 4 variables.
#' @seealso \code{\link{mcmcPottsNoData}}
"res5"

#' Results from fitting the surrogate model in Stan.
#' 
#' @format A \code{list} containing 3 variables:
#' \describe{
#'   \item{summary}{summary of the model fit for all chains merged}
#'   \item{c_summary}{summary for each HMC chain individually}
#'   \item{tm}{time taken to fit the model}
#' }
#' @seealso \code{\link[rstan]{summary,stanfit-method}}
"fitStan"

#' Results from fitting the hidden Potts model using
#' the approximate exchange algorithm (AEA).
#' 
#' @format A \code{list} containing 8 variables.
#' @seealso \code{\link{mcmcPotts}}
"resAEA"

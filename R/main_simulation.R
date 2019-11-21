####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions for simulating Gaussian Processes (GPs).
###
###   Description:
###
###
####################################################################################

#' @title Naive brute force simulation of gaussian processes
#'
#' @description Naive brute force simulation of gaussian processes
#'
#' @param n the number of locations
#' @param p the dimension of domain
#' @param loc.method the method to choose locations. The default method \code{"random"} uses randomly selected locations. The alternative \code{"grid"} uses grids as locations.
#' @param covmodel covariance model
#' @param covparms covariance parameters as a vector
#' @param nugget either a single (constant) nugget or a vector of nugget terms for the observations
#'
#' @return \code{simulate_gp_brute} returns
#' \itemize{
#'       \item{\code{n}: }{the number of locations}
#'       \item{\code{p}: }{the dimension of domain}
#'       \item{\code{locs}: }{a \code{n} x \code{p} matrix of locations. Each row of \code{locs} contains a location.}
#'       \item{\code{sigma}: }{covariance matrix with respect to \code{locs}}
#'       \item{\code{y}: }{a \code{n}-dimensional vector of the user-specific gaussian process generated at \code{locs}}
#' }
#' @export
#'
#' @examples
#' out <- simulate_gp_brute(n = 5^2, p = 2, loc.method = "grid", covmodel = diag(5^2), covparms = NULL, nugget = 0)
#' out
simulate_gp_brute <- function(n, p, loc.method = "random", covmodel, covparms = NULL, nugget = 0) {

  if(loc.method == "grid") {
    grid.oneside    <- seq(from = 0, to = 1, length.out = n^(1/p))
    locs            <- expand.grid( as.data.frame(matrix(grid.oneside, nrow = length(grid.oneside), ncol = p, byrow = F)) )

    if(n != nrow(locs)) {
      n             <- nrow(locs)
      message("Warning: The arguments n and p are not compatible. n will be ", n, ".")
    }

  } else if(loc.method == "random") {
    locs            <- as.data.frame(matrix(runif(n * p, 0, 1), n, p))
  } else {
    message("Warning: The argument loc.method is not valid. Locations are randomly selected.")
    locs            <- as.data.frame(matrix(runif(n * p, 0, 1), n, p))
  }

  if(is.function(covmodel)) {
    sigma           <- covmodel(locs1 = locs, locs2 = locs, covparms = covparms)
  } else if(is.matrix(covmodel)) {
    sigma           <- covmodel
  } else {
    message("Warning: The argument covmodel is not valid. The covariance matrix will be an identity matrix.")
    sigma           <- diag(n)
  }

  y                 <-  as.numeric(t(chol(sigma)) %*% rnorm(n))

  return(list(n = n, p = p, locs = locs, sigma = sigma, y = y))
}

#' Title
#'
#' @param temp1
#' @param temp2
#'
#' @return
#' @export
#'
#' @examples
simulate_gp <- function(temp1, temp2) {
  return(temp1 + temp2)
}

.simulate_gp_stationary <- function() {
  return(1)
}

.simulate_gp_nonstationary <- function() {
  return(1)
}

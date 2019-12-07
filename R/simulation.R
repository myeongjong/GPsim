####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions for simulating Gaussian Processes (GPs).
###
####################################################################################


#' @title Simulation of gaussian processes
#'
#' @description This function returns a realization of a user-specific gaussian process. This function is based on the fact that a gaussian process is fully specified by a mean function and a covariance matrix.
#'
#' @param locs If \code{locs} is a numeric matrix, then \code{locs} is a location matrix with \code{n} rows and \code{p} columns.
#'             If \code{locs} is \code{"grid"}, then a location matrix is generated using equidistant grids with the number \code{n} of locations and the dimension \code{p} of domain \code{[0,1]^p}.
#'             If \code{locs} is \code{"random"}, then a location matrix is randomly generated with the number \code{n} of locations and the dimension \code{p} of domain \code{[0,1]^p}.
#'             At \code{NULL} by default
#' @param n Number of locations. At \code{NULL} by default
#' @param p Dimension of domain \code{[0,1]^p}. At \code{NULL} by default
#' @param meanmodel If \code{meanmodel} is a function, then \code{meanmodel} is a mean function.
#'                  If \code{meanmodel} is a vector, then \code{meanmodel} is a mean vector.
#'                  By default, it is the zero function.
#' @param meanparms A numerical vector with mean parameters. It must be compatible with the argument \code{meanmodel}. At \code{NULL} by default
#' @param covmodel If \code{covmodel} is a function, then \code{covmodel} is a covariance function.
#'                 If \code{covmodel} is a matrix, then \code{covmodel} is a covariance matrix.
#'                 By default, it is the identity matrix.
#' @param covparms A numerical vector with covariance parameters. It must be compatible with the argument \code{covmodel}. At \code{NULL} by default
#' @param pivot Logical indicating if pivoting is to be used when factorizing a covariance matrix.  At \code{FALSE} by default
#' @param correction An argument specifying a correction method for the cholesky factorization of a covariance matrix. At \code{NULL} by default.
#'                   If correction is \code{NULL}, then the built-in R function \code{chol} is used.
#'                   If correction is \code{"qr"}, then the built-in R function \code{qr} is used.
#'                   If correction is \code{"diag"}, then \code{C + diag(tol, n)} is used instead of a covariance matrix \code{C}.
#'                   Other correction methods \code{"type-I"}, \code{"type-II"}, \code{"eigen-I"}, \code{"eigen-II"}, \code{"GMW81"}, and \code{"SE99"} are refered to Fang and O'leary (2008).
#' @param tol Numerical tolerance. At \code{.Machine$double.eps} by default
#' @param seed An integer specifying whether and how the random number generator should be initialized. The default, \code{NULL} will not change the random generator state.
#'
#' @return \code{simulate_gp} returns
#'     \itemize{
#'         \item{\code{seed}: } The random seed
#'         \item{\code{n}: } The number of locations
#'         \item{\code{p}: } The dimension of domain \code{[0,1]^p}
#'         \item{\code{locs}: } The location matrix with \code{n} rows and \code{p} columns. Each row of \code{locs} gives a point in \code{[0, 1]^p}
#'         \item{\code{meanvec}: } The mean vector induced by \code{meanmodel}
#'         \item{\code{covmat}: } The covariance matrix induced by \code{covmodel}
#'         \item{\code{y}: } A \code{n}-dimensional vector of the user-specific gaussian process generated at \code{locs}
#'     }
#'
#' @references Fang, Haw-ren, and Dianne P. O’leary. "Modified Cholesky algorithms: a catalog with new approaches." Mathematical Programming 115.2 (2008): 319-349.
#'
#' @export
#'
#' @examples
#' out <- simulate_gp(locs = NULL, n = 5^2, p = 2)
#' out$y # a realization of the gaussian process with a mean rep(0, n) and a covariance diag(1, n)
simulate_gp <- function(locs = NULL, n = NULL, p = NULL, meanmodel = function(locs, meanparms) rep(0, nrow(locs)), meanparms = NULL, covmodel = function(locs, covparms) diag(1, nrow(locs)), covparms = NULL, pivot = FALSE, correction = NULL, tol = .Machine$double.eps, seed = NULL) {

  # set a random seed
  if(!is.null(seed)) {
    seed.old <- .Random.seed
    on.exit( { .Random.seed <<- seed.old } )
    set.seed(seed)
  }

  # check arguments: locs, n, and p
  locsnp    <- .checkargs_locsnp(locs = locs, n = n, p = p)

  # check arguments: meanmodel and meanparms
  meanlocs  <- .checkargs_meanmodel(meanmodel = meanmodel, meanparms = meanparms, locs = locsnp$locs, n = locsnp$n)

  # check arguments: covmodel and covparms
  sigma     <- .checkargs_covmodel(covmodel = covmodel, covparms = covparms, locs = locsnp$locs, n = locsnp$n)

  # generate a realization of the gaussian process
  gp        <- .process_gp(meanlocs = meanlocs, covlocs = sigma, n = locsnp$n, pivot = pivot, correction = correction, tol = tol)

  # return
  return(list(seed = seed, n = locsnp$n, p = locsnp$p, locs = locsnp$locs, meanvec = meanlocs, covmat = sigma, y = gp$y))
}

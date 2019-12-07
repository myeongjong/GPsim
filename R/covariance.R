####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes two R functions for implementing covariance functions: one is an isotropic exponential covariance function and the other is an anisotropic exponential covariance function.
###
####################################################################################


#' @title Isotropic exponential covariance function
#'
#' @description From a location matrix (\code{locs}) and a vector with covariance parameters (\code{covparms}), this function returns an isotropic exponential covariance matrix which is one of the simplest covariance matrices.
#'
#' @param locs A matrix with \code{n} rows and \code{p} columns. Each row of locs gives a point in R^p
#' @param covparms A vector with covariance parameters in the form (variance, range)
#'
#' @section Parametrization: The covariance parameter vector is (variance, range) = \eqn{(\sigma^2 , r)}. The form of the covariance is \deqn{ C(x, y) = \sigma^2 exp( || x - y || / r)} where \eqn{x} and \eqn{y} are locations in R^p.
#'
#' @return A matrix with \code{n} rows and \code{n} columns, with the (i, j) entry containing the isotropic exponenital covariance between observations \code{locs[i, ]} and \code{locs[j, ]}
#' @export
#'
#' @examples
#' # grid locations
#' cov_expo_iso(locs = expand.grid(c(0.25, 0.75), c(0.25, 0.75)), covparms = c(1, 0.1))
#'
#' # randomly selected locations
#' cov_expo_iso(locs = matrix(runif(8), 4, 2), covparms = c(1, 0.1))
cov_expo_iso <- function(locs, covparms) {

  covparms[1] * exp(-fields::rdist(locs) / covparms[2])
}


#' @title Anisotropic exponential covariance function
#'
#' @description From a location matrix (\code{locs}) and a vector with covariance parameters (\code{covparms}), this function returns an anisotropic exponential covariance matrix.
#'
#' @param locs A matrix with \code{n} rows and \code{p} columns. Each row of locs gives a point in R^d
#' @param covparms A vector with covariance parameters in the form (variance, range, degree of anisotropy)
#'
#' @section Parametrization: The covariance parameter vector is (variance, range, degree of anisotropy) = \eqn{(\sigma^2 , r, \alpha)}. The form of the covariance is \deqn{ C(x, y) = \sigma^2 exp( || A ( x - y ) || / r)} where \eqn{x} and \eqn{y} are locations in R^p and A is a diagonal matrix \eqn{diag( \sqrt \alpha , 1 , ... , 1 )}.
#'
#' @return A matrix with \code{n} rows and \code{n} columns, with the (i, j) entry containing the anisotropic exponenital covariance between observations \code{locs[i, ]} and \code{locs[j, ]}
#' @export
#'
#' @examples
#' # grid locations
#' cov_expo_aniso(locs = expand.grid(c(0.25, 0.75), c(0.25, 0.75)), covparms = c(1, 0.1, 10))
#'
#' # randomly selected locations
#' cov_expo_aniso(locs = matrix(runif(8), 4, 2), covparms = c(1, 0.1, 10))
cov_expo_aniso <- function(locs, covparms) {

  covparms[1] * exp(-fields::rdist(cbind(locs[ ,1] * covparms[3], locs[ ,-1])) / covparms[2])
}

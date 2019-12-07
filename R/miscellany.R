####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions for
###
###   Description:
###
###
####################################################################################


#' @title Kullback Leiber (KL) divergence between two multivariate normal distributions
#'
#' @description From two mean vectors (\code{mu0}, \code{mu1}) and two covariance matrices (\code{covmat0}, \code{covmat1}), this function returns KL divergence between \code{N(mu0,covmat0)} and \code{N(mu1,covmat1)}.
#'
#' @param covmat0 A covariance matrix of the first multivariate normal distribution
#' @param covmat1 A covariance matrix of the second multivariate normal distribution
#' @param mu0 A mean vector of the first multivariate normal distribution
#' @param mu1 A mean vector of the second multivariate normal distribution
#'
#' @return KL divergence between \code{N(mu0,covmat0)} and \code{N(mu1,covmat1)}
#' @export
#'
#' @examples
#' covmat0 <- diag(3)
#' covmat1 <- matrix(c(1, 0.9, 0.7, 0.9, 1, 0.4, 0.7, 0.4, 1), 3, 3)
#'
#' kldiv(covmat0 = covmat0, covmat1 = covmat1)
kldiv <- function(covmat0, covmat1, mu0 = rep(0, nrow(covmat0)), mu1 = rep(0, nrow(covmat0))) {

  # function to compute KL divergence
  # multivariate normal KL divergence between true distribution N(mu0,covmat0) and approx distr N(mu1,covmat1)

  if(!is.matrix(covmat0) | !is.matrix(covmat1)) stop("Both covmat0 and covmat1 must be matrices.")

  n             <- nrow(covmat0)
  chol0         <- t(chol(covmat0))
  chol1         <- t(chol(covmat1))

  # trace term
  M             <- solve(chol1, chol0)
  traceterm     <- sum(M^2) - n

  # log det term
  logdetterm    <- 2 * sum( log(diag(chol1)) ) - 2 * sum( log(diag(chol0)) )

  # mean term (equal to zero if means are zero)
  meandiff      <- mu1 - mu0
  # t(meandiff) %*% solve(covmat1,meandiff)
  temp          <- forwardsolve(chol1, meandiff)
  meanterm      <- sum(temp^2)

  # kldiv
  kldiv         <- 1/2 * (traceterm + logdetterm + meanterm)

  return(kldiv)

}

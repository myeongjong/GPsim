####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions for simulating Gaussian Processes (GPs).
###
###   Description:
###
###
####################################################################################


#' @title Brute-force simulation of gaussian processes
#'
#' @description Brute-force simulation of gaussian processes
#'
#' @param locs locations matrix or how to select locations from the unit box
#' @param n number of locations
#' @param p dimension of the unit box
#' @param meanmodel mean model
#' @param meanparms mean parameters as a vector
#' @param covmodel covariance model
#' @param covparms covariance parameters as a vector
#' @param seed random seed
#'
#' @return \code{simulate_gp_brute} returns a list containing
#' \itemize{
#'       \item{seed: }{random seed}
#'       \item{n: }{number of locations}
#'       \item{p: }{dimension of the unit box}
#'       \item{locs: }{a \code{n} x \code{p} matrix of locations}
#'       \item{meanvec: }{mean values at locations}
#'       \item{covmat: }{covariance matrix}
#'       \item{y: }{generated gaussian process at locations}
#' }
#' @export
#'
#' @examples
#' # The following arguments give the exact same result:
#'
#' simulate_gp_brute(locs = 'random', n = 3, p = 2, meanmodel = function(loc, meanparms) 0, meanparms = NULL, covmodel = function(loc1, loc2, covparms) ifelse(identical(loc1, loc2), 1, 0), covparms = NULL, seed = 1)$y
#'
#' simulate_gp_brute(locs = 'random', n = 3, p = 2, meanmodel = matrix(0, 3, 1), meanparms = NULL, covmodel = diag(3), covparms = NULL, seed = 1)$y
#'
#' simulate_gp_brute(locs = 'random', n = 3, p = 2, meanmodel = matrix(0, 1, 3), meanparms = NULL, covmodel = diag(3), covparms = NULL, seed = 1)$y
#'
#' simulate_gp_brute(locs = 'random', n = 3, p = 2, meanmodel = as.data.frame(matrix(0, 3, 1)), meanparms = NULL, covmodel = as.data.frame(diag(3)), covparms = NULL, seed = 1)$y
#'
#' simulate_gp_brute(locs = 'random', n = 3, p = 2, meanmodel = as.data.frame(matrix(0, 1, 3)), meanparms = NULL, covmodel = as.data.frame(diag(3)), covparms = NULL, seed = 1)$y
#'
#' simulate_gp_brute(locs = 'random', n = 3, p = 2, meanmodel = matrix(0, 3, 1), meanparms = NULL, covmodel = rep(1, 3), covparms = NULL, seed = 1)$y
#'
#' simulate_gp_brute(locs = 'random', n = 3, p = 2, meanmodel = matrix(0, 1, 3), meanparms = NULL, covmodel = 'identity', covparms = NULL, seed = 1)$y
simulate_gp_brute <- function(locs = NULL, n = NULL, p = NULL, meanmodel = function(loc, meanparms) 0, meanparms = NULL, covmodel = function(loc1, loc2, covparms) ifelse(identical(loc1, loc2), 1, 0), covparms = NULL, seed = NULL) {

  if(!is.null(seed)) {
    seed.old <- .Random.seed
    on.exit( { .Random.seed <<- seed.old } )
    set.seed(seed)
  }

  locsnp <- .input_check_locs_n_p(locs, n, p)

  if(is.function(meanmodel)) {

    meanlocs <- apply(X = locsnp$locs, MARGIN = 1, FUN = function(x) meanmodel(loc = x, meanparms = meanparms))

  } else if(is.matrix(meanmodel) | is.data.frame(meanmodel)) {

    if(identical(as.numeric(dim(meanmodel)), c(locsnp$n, 1)) | identical(as.numeric(dim(meanmodel)), c(1, locsnp$n))) {
      meanmodel <- as.matrix(meanmodel)
      meanlocs <- as.numeric(meanmodel)
    } else {
      stop("The mean model (meanmodel) is not compatible with the locations (locs).")
    }

  } else if(is.numeric(meanmodel) | is.integer(meanmodel)) {

    if(length(meanmodel) == locsnp$n) {
      meanlocs <- as.numeric(meanmodel)
    } else {
      stop("The mean model (meanmodel) is not compatible with the locations (locs).")
    }

  } else {

    stop("The mean model (meanmodel) is not valid. Please refer the description.")

  }

  if(is.function(covmodel)) {

    ind <- expand.grid(seq(locsnp$n), seq(locsnp$n))
    covlocs <- apply(X = ind, MARGIN = 1, FUN = function(x) covmodel(loc1 = as.numeric(locsnp$locs[x[1], ]), loc2 = as.numeric(locsnp$locs[x[2], ]), covparms = covparms))
    covlocs <- matrix(covlocs, locsnp$n, locsnp$n, byrow = F)

  } else if(is.matrix(covmodel) | is.data.frame(covmodel)) {

    if(identical(as.numeric(dim(covmodel)), c(locsnp$n, locsnp$n))) {
      covlocs <- as.matrix(covmodel)
    } else {
      stop("The covariance model (covmodel) is not compatible with the locations (locs).")
    }

  } else if(is.numeric(covmodel) | is.integer(covmodel)) {

    if(length(covmodel) == locsnp$n) {
      covlocs <- diag(as.numeric(covmodel))
    } else {
      stop("The covariance model (covmodel) is not valid. Please refer the description.")
    }

  } else if(covmodel == "identity"){

    covlocs <- diag(locsnp$n)

  } else {

    stop("The covariance model (covmodel) is not valid. Please refer the description.")

  }

  y <- meanlocs + as.numeric(t(chol(covlocs)) %*% rnorm(locsnp$n))

  return(list(seed = seed, n = locsnp$n, p = locsnp$p, locs = locsnp$locs, meanvec = meanlocs, covmat = covlocs, y = y))
}

simulate_gp <- function(temp1, temp2) {
  return(temp1 + temp2)
}


.simulate_gp_stationary <- function() {
  return(1)
}

.simulate_gp_nonstationary <- function() {
  return(1)
}

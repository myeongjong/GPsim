####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions for simulating Gaussian Processes (GPs).
###
###   Description:
###
###
####################################################################################


simulate_gp_brute <- function(locs = NULL, n = NULL, p = NULL, meanmodel = function(loc, meanparms) 0, meanparms = NULL, covmodel = function(loc1, loc2, covparms) ifelse(identical(loc1, loc2), 1, 0), covparms = NULL, seed = NULL, correction = NULL, tol = 1e8 * .Machine$double.eps) {

  if(!is.null(seed)) {
    seed.old <- .Random.seed
    on.exit( { .Random.seed <<- seed.old } )
    set.seed(seed)
  }

  locsnp    <- .checkargs_locsnp(locs = locs, n = n, p = p)

  meanlocs  <- .checkargs_meanmodel(meanmodel = meanmodel, meanparms = meanparms, locs = locsnp$locs, n = locsnp$n, p = locsnp$p)

  sigma     <- .checkargs_covmodel(covmodel = covmodel, covparms = covparms, locs = locsnp$locs, n = locsnp$n, p = locsnp$p)

  gp        <- .process_gp(meanlocs = meanlocs, covlocs = sigma, n = locsnp$n, correction = correction, tol = tol)

  return(list(seed = seed, n = locsnp$n, p = locsnp$p, locs = locsnp$locs, meanvec = meanlocs, covmat = sigma, y = gp$y))
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

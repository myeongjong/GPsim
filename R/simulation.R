####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions for simulating Gaussian Processes (GPs).
###
###   Description:
###
###
####################################################################################


simulate_gp <- function(locs = NULL, n = NULL, p = NULL, meanmodel = function(locs, meanparms) rep(0, nrow(locs)), meanparms = NULL, covmodel = function(locs, covparms) diag(1, nrow(locs)), covparms = NULL, pivot = FALSE, correction = NULL, tol = .Machine$double.eps, seed = NULL) {

  if(!is.null(seed)) {
    seed.old <- .Random.seed
    on.exit( { .Random.seed <<- seed.old } )
    set.seed(seed)
  }

  locsnp    <- .checkargs_locsnp(locs = locs, n = n, p = p)

  meanlocs  <- .checkargs_meanmodel(meanmodel = meanmodel, meanparms = meanparms, locs = locsnp$locs, n = locsnp$n, p = locsnp$p)

  sigma     <- .checkargs_covmodel(covmodel = covmodel, covparms = covparms, locs = locsnp$locs, n = locsnp$n, p = locsnp$p)

  gp        <- .process_gp(meanlocs = meanlocs, covlocs = sigma, n = locsnp$n, pivot = pivot, correction = correction, tol = tol)

  return(list(seed = seed, n = locsnp$n, p = locsnp$p, locs = locsnp$locs, meanvec = meanlocs, covmat = sigma, y = gp$y))
}

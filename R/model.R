####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions for implementing baseline models.
###
###   Description:
###
###
####################################################################################


corrvecchia_knownCovparms <- function(locs, m, ordering = "maxmin", coordinate = NULL, def.dist = NULL, ordering.method = "correlation", initial.pt = NULL, conditioning = "NN", conditioning.method = "correlation", covmodel, covparms) {

  locs  <- as.matrix(locs)
  p     <- ncol(locs)
  n     <- nrow(locs)

  if(ordering.method == "correlation" | conditioning.method == "correlation") {
    rho         <- .correlation(locs = locs, covmodel = covmodel, covparms = covparms, def.dist = def.dist)
  }

  # ordering
  if(ordering == "coord") {
    if(is.null(coordinate)) coordinate <- seq(p)
    ord         <- .order_coordinate(locs = locs, coordinate = coordinate)
  } else if(ordering == "maxmin" & ordering.method == "euclidean") {
    ord         <- .order_maxmin_euclidean(locs = locs)
  } else if(ordering == "maxmin" & ordering.method == "correlation") {
    ord         <- .order_maxmin_correlation(locs = locs, dinv = rho, covmodel = covmodel, covparms = covparms, initial.pt = initial.pt)
  } else {
    stop("Please check the ordering method.")
  }

  locsord       <- locs[ord, , drop = FALSE]

  if(is.matrix(covmodel)) {
    covmodel    <- covmodel[ord, ord]
  }

  # conditioning
  if(conditioning.method == "euclidean") {
    cond.sets   <- GpGp::find_ordered_nn(locs = locsord, m = m)
  } else if(conditioning.method == "correlation") {
    rho         <- rho[ord, ord]
    cond.sets   <- .conditioning_nn(m = m, dist.matrix = 1 - rho)
  } else {
    stop("Please check the conditioning method.")
  }

  Cond          <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
  obs           <- rep(TRUE, n)
  U.prep        <- .U_sparsity(locsord, cond.sets, obs, Cond)

  # output
  vecchia.approx <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord, ord.pred='general', U.prep = U.prep, cond.yz = 'false', ordering = ordering, ordering.method = ordering.method, conditioning = 'NN', conditioning.method = conditioning.method)
  return(vecchia.approx)
}

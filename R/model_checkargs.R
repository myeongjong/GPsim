####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions for
###
###   Description:
###
###
####################################################################################


.checkargs_locsm_model <- function(locs, m) {

  # the argument: locs
  if( !(is.matrix(locs) | is.data.frame(locs)) ) stop("Locations (locs) must be either a matrix or data frame.")

  locs <- as.matrix(locs)

  if(!is.numeric(locs)) stop("Locations (locs) must be numeric.")

  # the argument: m
  if(!is.numeric(m)) stop("The size of conditioning sets (m) is not numeric.")

  if(length(m) > 1) stop("The size of conditioning sets (m) is not a single number.")


  # compatibility
  if(nrow(locs) <= m) stop("The size of conditioning sets (m) must be less than the number of locations.")

  return(0)
}


.checkargs_ordering_model <- function(ordering, coordinate, dist.ordering, p) {

  if(ordering == "coord" & dist.ordering == "euclidean") {

    if(is.null(coordinate)) {
      warning("Coordinate must be specified when using coordinate-based ordering. For convenience, every coordinates are used.")
      coordinate <- seq(p)
    }

    if(!is.numeric(coordinate)) stop("The argument coordinate must be numeric.")

    if(length(coordinate) != length(unique(coordinate))) stop("The entries of the argument coordinate must be distinct.")

    if(!(min(coordinate) >= 1 & max(coordinate) <= p)) stop("The entries of the argument coordinate must be in between 1 and the dimension of domain.")

  } else if(ordering == "coord" & dist.ordering == "correlation") {

    stop("Coordinate-based ordering (ordering = coord) is not compatible with the correlation-based distance. Please select maxmin ordering (ordering = maxmin).")

  } else if(ordering == "maxmin" & dist.ordering == "euclidean") {

    if(!is.null(coordinate)) warning("The argument coordinate is not used.")

  } else if(ordering == "maxmin" & dist.ordering == "correlation") {

    if(!is.null(coordinate)) warning("The argument coordinate is not used.")

  } else {

    stop("Arguments related to ordering is invalid. Please check the arguments: ordering, coordinate, dist.ordering.")

  }

  return(0)
}


.checkargs_conditioning_model <- function(dist.conditioning) {

  if(dist.conditioning == "euclidean" | dist.conditioning == "correlation") {

    # message("Note: This function only consider nearest neighbor conditioning sets so far.")

  } else {

    stop("An argument related to conditioning is invalid. Please check the argument dist.conditioning.")

  }

  return(0)
}


.checkargs_covmodel_model <- function(covmodel, covparms, locs, n, ncheck = 10) {

  if(is.function(covmodel)) {

    ncheck <- min(n, ncheck)
    covlocs <- tryCatch(covmodel(locs[1:ncheck, ], covparms), error = function(e) "Oops")

    if(identical(covlocs, "Oops")) stop("Either the covariance function (covmodel) or its parameters (covparms) do not work properly.")

  } else if(is.matrix(covmodel) | is.data.frame(covmodel)) {

    if( !(nrow(covmodel == n) & ncol(covmodel) == n) ) stop("The covariance matrix (covmodel) is not compatible with the locations (locs).")

  } else {

    stop("The covariance model (covmodel) is not valid. Please refer the description.")

  }

  return(0)
}

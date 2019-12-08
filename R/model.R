####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions for implementing baseline models such as the euclidean based and the correlation based Vecchia approximation.
###
####################################################################################


#' @title Vecchia approximation with Euclidean and correlation-based distances
#'
#' @description The correlation-based Vecchia approximation is nothing but the Vecchia approximation with a correlation-based distance.
#'              It is equivalent to the Vecchia approximation with Euclidean distance for isotropic covariance function cases which are popular in application.
#'              If offers an automatic strategy even when Euclidean distance is not applicable (e.g. text data).
#'
#' @param locs A matrix with \code{n} rows and \code{p} columns. Each row of locs gives a point in \code{[0, 1]^p}
#' @param m Number of nearby points to condition on (the size of conditioning sets)
#' @param ordering "coord" or "maxmin."
#'                 If \code{ordering} is "coord," then coordinate-based ordering method is used to order the locations.
#'                 If \code{ordering} is "maxmin," then maxmin ordering method is used to order the locations.
#' @param coordinate a numeric vector of coordinates
#' @param dist.ordering "euclidean" or "correlation."
#'                      If \code{dist.ordering} is "euclidean," then euclidean distance is used to order the locations.
#'                      If \code{dist.ordering} is "correlation," then correlation based distance 1-rho is used to order the locations.
#' @param dist.conditioning "euclidean" or "correlation."
#'                          If \code{dist.conditioning} is "euclidean," then euclidean distance is used to construct conditioning sets.
#'                          If \code{dist.conditioning} is "correlation," then correlation based distance 1-rho is used to construct conditioning sets.
#' @param covmodel If \code{covmodel} is a function, then \code{covmodel} is a covariance function.
#'                 If \code{covmodel} is a matrix, then \code{covmodel} is a covariance matrix.
#'                 Please use \code{covparms = c(1)} if \code{covmodel} is a correlation matrix.
#' @param covparms A numerical vector with covariance parameters. It must be compatible with the argument \code{covmodel}. At \code{NULL} by default
#'
#' @return An object that specifies the Vecchia approximation for later use in likelihood evaluation or prediction. We are doing research on this.
#'
#' @references Katzfuss, Matthias, and Joseph Guinness. "A general framework for Vecchia approximations of Gaussian processes." arXiv preprint arXiv:1708.06302 (2017).
#'
#' @export
#'
#' @examples
#' n    <- 15^2
#' m    <- 10
#' locs <- matrix(runif(n * 2, 0, 1), n, 2)
#'
#' covparms   <- c(1, 0.1, 10)
#' sigma      <- cov_expo_aniso(locs = locs, covparms = covparms)
#'
#' out.euclidean <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", coordinate = NULL, dist.ordering = "euclidean", dist.conditioning = "euclidean", covmodel = cov.aniso, covparms = covparms)
#' out.correlation <- corrvecchia_knownCovparms(locs = locs, m = m, ordering = "maxmin", coordinate = NULL, dist.ordering = "correlation", dist.conditioning = "correlation", covmodel = cov.aniso, covparms = covparms)
#'
#' out.euclidean$ord
#' out.correlation$ord
corrvecchia_knownCovparms <- function(locs, m, ordering = "maxmin", coordinate = NULL, dist.ordering = "correlation", dist.conditioning = "correlation", covmodel, covparms = NULL) {

  .checkargs_locsm_model(locs = locs, m = m)

  locs  <- as.matrix(locs)
  n     <- nrow(locs)
  p     <- ncol(locs)

  .checkargs_ordering_model(ordering = ordering, coordinate = coordinate, dist.ordering = dist.ordering, p = p)

  .checkargs_conditioning_model(dist.conditioning = dist.conditioning)

  .checkargs_covmodel_model(covmodel = covmodel, covparms = covparms, locs = locs, n = n, ncheck = 10)

  # calculate a correlation matrix
  if(dist.ordering == "correlation" | dist.conditioning == "correlation") {
    rho         <- .correlation(locs = locs, covmodel = covmodel, covparms = covparms, def.dist = NULL)
  }

  # ordering
  if(ordering == "coord") {

    if(is.null(coordinate)) coordinate <- seq(p)
    ord         <- .order_coordinate(locs = locs, coordinate = coordinate)

  } else if(ordering == "maxmin" & dist.ordering == "euclidean") {

    ord         <- .order_maxmin_euclidean(locs = locs)

  } else if(ordering == "maxmin" & dist.ordering == "correlation") {

    ord         <- .order_maxmin_correlation(locs = locs, dinv = rho, covmodel = covmodel, covparms = covparms, initial.pt = NULL)

  } else {

    stop("Please check the ordering method.")

  }

  locsord       <- locs[ord, , drop = FALSE]
  if(is.matrix(covmodel)) covmodel <- covmodel[ord, ord]

  # conditioning
  if(dist.conditioning == "euclidean") {

    cond.sets   <- GpGp::find_ordered_nn(locs = locsord, m = m)

  } else if(dist.conditioning == "correlation") {

    rho         <- rho[ord, ord]
    cond.sets   <- .conditioning_nn(m = m, dist.matrix = 1 - rho)

  } else {

    stop("Please check the conditioning method.")

  }

  cond          <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); cond[!is.na(cond.sets)] <- TRUE
  obs           <- rep(TRUE, n)
  U.prep        <- .U_sparsity(locsord, cond.sets, obs, cond)

  # output
  vecchia.approx <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord, ord.pred='general', U.prep = U.prep, cond.yz = 'false', ordering = ordering, dist.ordering = dist.ordering, conditioning = 'NN', dist.conditioning = dist.conditioning)
  return(vecchia.approx)
}


# corrvecchia_knownCovparms_old <- function(locs, m, ordering = "maxmin", coordinate = NULL, def.dist = NULL, ordering.method = "correlation", initial.pt = NULL, conditioning.method = "correlation", covmodel, covparms) {
#
#   locs  <- as.matrix(locs)
#   p     <- ncol(locs)
#   n     <- nrow(locs)
#
#   if(ordering.method == "correlation" | conditioning.method == "correlation") {
#     rho         <- .correlation(locs = locs, covmodel = covmodel, covparms = covparms, def.dist = def.dist)
#   }
#
#   # ordering
#   if(ordering == "coord") {
#     if(is.null(coordinate)) coordinate <- seq(p)
#     ord         <- .order_coordinate(locs = locs, coordinate = coordinate)
#   } else if(ordering == "maxmin" & ordering.method == "euclidean") {
#     ord         <- .order_maxmin_euclidean(locs = locs)
#   } else if(ordering == "maxmin" & ordering.method == "correlation") {
#     ord         <- .order_maxmin_correlation(locs = locs, dinv = rho, covmodel = covmodel, covparms = covparms, initial.pt = initial.pt)
#   } else {
#     stop("Please check the ordering method.")
#   }
#
#   locsord       <- locs[ord, , drop = FALSE]
#
#   if(is.matrix(covmodel)) {
#     covmodel    <- covmodel[ord, ord]
#   }
#
#   # conditioning
#   if(conditioning.method == "euclidean") {
#     cond.sets   <- GpGp::find_ordered_nn(locs = locsord, m = m)
#   } else if(conditioning.method == "correlation") {
#     rho         <- rho[ord, ord]
#     cond.sets   <- .conditioning_nn(m = m, dist.matrix = 1 - rho)
#   } else {
#     stop("Please check the conditioning method.")
#   }
#
#   Cond          <- matrix(NA, nrow(cond.sets), ncol(cond.sets)); Cond[!is.na(cond.sets)] <- TRUE
#   obs           <- rep(TRUE, n)
#   U.prep        <- .U_sparsity(locsord, cond.sets, obs, Cond)
#
#   # output
#   vecchia.approx <- list(locsord = locsord, obs = obs, ord = ord, ord.z = ord, ord.pred='general', U.prep = U.prep, cond.yz = 'false', ordering = ordering, ordering.method = ordering.method, conditioning = 'NN', conditioning.method = conditioning.method)
#   return(vecchia.approx)
# }

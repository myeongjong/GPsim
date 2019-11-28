####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions for
###
###   Description:
###
###
####################################################################################


.process_gp <- function(meanlocs, covlocs, n, correction = NULL, tol = 1e8 * .Machine$double.eps) {

  if(is.null(correction)) {

    covfactor             <- tryCatch(chol(covlocs), error = function(e) "Oops")

    if(identical(covfactor, "Oops")) stop("Cholesky decomposition cannnot be applied to the covariance matrix. Please specify the correction method (correction).")

  } else if(correction == "qr") {

    eigendecomp           <- eigen(covlocs)
    covfactor             <- tryCatch(qr.R(qr( diag(sqrt(eigendecomp$values)) %*% t(eigendecomp$vectors) )), error = function(e) "Oops")

    if(identical(covfactor, "Oops")) stop("QR decomposition cannnot be applied. Please use another correction method (correction).")

  } else if(correction == "diagonal") {

    covfactor             <- tryCatch(chol(covlocs + diag(tol, n)), error = function(e) "Oops")

    if(identical(covfactor, "Oops")) stop("Cholesky decomposition cannnot be applied to the modified covariance matrix. Please use larger tolerance value (tol).")

  } else if(correction == "type1") {

    eigendecomp           <- eigen(covlocs)
    eigendecomp$values    <- ifelse(abs(eigendecomp$values) > tol, abs(eigendecomp$values), tol)
    covfactor             <- tryCatch(qr.R(qr( diag(sqrt(eigendecomp$values)) %*% t(eigendecomp$vectors) )), error = function(e) "Oops")

    if(identical(covfactor, "Oops")) stop("QR decomposition cannnot be applied. Please use larger tolerance value (tol).")

  } else if(correction == "type2") {

    eigendecomp           <- eigen(covlocs)
    eigendecomp$values    <- ifelse(eigendecomp$values > tol, eigendecomp$values, tol)
    covfactor             <- tryCatch(qr.R(qr( diag(sqrt(eigendecomp$values)) %*% t(eigendecomp$vectors) )), error = function(e) "Oops")

    if(identical(covfactor, "Oops")) stop("QR decomposition cannnot be applied. Please use larger tolerance value (tol).")

  } else {

    stop("Not yet.")

  }

  covlocs.modified <- t(covfactor) %*% covfactor

  err.decomp <- sqrt(sum( (covlocs - covlocs.modified)^2 ))
  message("The decomposition error of the covariance matrix is ", err.decomp, " with respect to the Frobenius norm.")

  y <- meanlocs + as.numeric(t(covfactor) %*% rnorm(n))

  return(list(y = y, meanlocs = meanlocs, covlocs = covlocs, covlocs.modified = covlocs.modified, correction = correction, tol = tol))
}

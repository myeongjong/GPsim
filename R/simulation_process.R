####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions for
###
###   Description:
###
###
####################################################################################


.process_gp <- function(meanlocs, covlocs, n, correction = NULL, tol = .Machine$double.eps) {

  if(is.null(correction)) {

    covfactor             <- tryCatch(chol(covlocs), error = function(e) "Oops")

    if(identical(covfactor, "Oops")) stop("Cholesky decomposition cannnot be applied to the covariance matrix. Please specify the correction method (correction).")

  } else if(correction == "qr") {

    eigendecomp           <- eigen(covlocs)
    covfactor             <- tryCatch(qr.R(qr( diag(sqrt(eigendecomp$values)) %*% t(eigendecomp$vectors), tol = tol )), error = function(e) "Oops")

    if(identical(covfactor, "Oops")) stop("QR decomposition cannnot be applied. Please use another correction method (correction).")

  } else if(correction == "diagonal") {

    covfactor             <- tryCatch(chol(covlocs + diag(tol, n)), error = function(e) "Oops")

    if(identical(covfactor, "Oops")) stop("Cholesky decomposition cannnot be applied to the modified covariance matrix. Please use larger tolerance value (tol).")

  } else if(correction == "type-I") {

    eigendecomp           <- eigen(covlocs)
    eigendecomp$values    <- ifelse(abs(eigendecomp$values) > tol, abs(eigendecomp$values), tol)
    covfactor             <- tryCatch(qr.R(qr( diag(sqrt(eigendecomp$values)) %*% t(eigendecomp$vectors) )), error = function(e) "Oops")

    if(identical(covfactor, "Oops")) stop("QR decomposition cannnot be applied. Please use larger tolerance value (tol).")

  } else if(correction == "type-II") {

    eigendecomp           <- eigen(covlocs)
    eigendecomp$values    <- ifelse(eigendecomp$values > tol, eigendecomp$values, tol)
    covfactor             <- tryCatch(qr.R(qr( diag(sqrt(eigendecomp$values)) %*% t(eigendecomp$vectors) )), error = function(e) "Oops")

    if(identical(covfactor, "Oops")) stop("QR decomposition cannnot be applied. Please use larger tolerance value (tol).")

  } else if(correction == "GMW81") {

    stop("Not yet.")

  } else if(correction == "SE90") {

    stop("Not yet.")

  } else if(correction == "SE99") {

    stop("Not yet.")

  } else if(correction == "GMW-I") {

    stop("Not yet.")

  } else if(correction == "GMW-II") {

    stop("Not yet.")

  } else if(correction == "SE-I") {

    stop("Not yet.")

  } else {

    stop("Not yet.")

  }

  covlocs.modified <- t(covfactor) %*% covfactor

  err.decomp <- sqrt(sum( (covlocs - covlocs.modified)^2 ))
  message("The decomposition error of the covariance matrix is ", err.decomp, " with respect to the Frobenius norm.")

  y <- meanlocs + as.numeric(t(covfactor) %*% rnorm(n))

  return(list(y = y, meanlocs = meanlocs, covlocs = covlocs, covlocs.modified = covlocs.modified, correction = correction, tol = tol))
}

.algorithm_SE <- function(A, delta) {

  n <- nrow(A)
  for(k in 1:(n-1)) {
    A[k, k] <- sqrt(A[k, k] + delta[k])
    A[(k+1):n, k] <- A[(k+1):n, k] / A[k, k]
    A[(k+1):n, (k+1):n] <- A[(k+1):n, (k+1):n] - tcrossprod(A[(k+1):n, k])
  }

  A[n, n] <- sqrt(A[n, n] + delta[n])
  A <- A * !upper.tri(A)
  return(A)
}

# B <- matrix(c(2, -1, 0, -1, 2, -1, 0, -1, 2), 3, 3)
# delta <- c(3, 1, 2)
# temp <- .algorithm_SE(B, delta)
# temp %*% t(temp)
#
# A <- matrix(c(1, 0.9, 0.7, 0.9, 1, 0.4, 0.7, 0.4, 1), 3, 3)
# delta <- c(3, 1, 2)
# temp <- .algorithm_SE(A, delta)
# temp %*% t(temp)

.algorithm_SE_LDL <- function(A, delta) {

  n <- nrow(A)

  L <- diag(n)
  D <- rep(NA, n)

  A.k <- A
  for(k in 1:(n-1)) {
    D[k]            <- A.k[1, 1] + delta[k]
    L[(k+1):n, k]   <- A.k[2:nrow(A.k), 1] / D[k]
    A.k             <- A.k[2:nrow(A.k), 2:nrow(A.k)] - tcrossprod(A.k[2:nrow(A.k), 1]) / D[k]
  }

  D[n] <- A.k[nrow(A.k), nrow(A.k)] + delta[n]
  return(list(Lmat = L, Dvec = D))
}

# B <- matrix(c(2, -1, 0, -1, 2, -1, 0, -1, 2), 3, 3)
# delta <- c(3, 1, 2)
# temp <- .algorithm_SE_LDL(B, delta)
# temp$Lmat %*% diag(temp$Dvec) %*% t(temp$Lmat)
#
# A <- matrix(c(1, 0.9, 0.7, 0.9, 1, 0.4, 0.7, 0.4, 1), 3, 3)
# delta <- c(3, 1, 2)
# temp <- .algorithm_SE_LDL(A, delta)
# temp$Lmat %*% diag(temp$Dvec) %*% t(temp$Lmat)

.algorithm_GMW81 <- function(A, beta, tol) {

  eigendecomp <- eigen(A)


}





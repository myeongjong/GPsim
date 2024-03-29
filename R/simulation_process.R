####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several small R functions which compose the simulate_gp() function.
###
####################################################################################


.process_gp <- function(meanlocs, covlocs, n, pivot = FALSE, correction = NULL, tol = .Machine$double.eps) {

  if(is.null(correction)) {

    covfactor             <- tryCatch(chol(x = covlocs, pivot = pivot, tol = tol), error = function(e) "Oops")
    if(identical(covfactor, "Oops")) stop("Cholesky decomposition cannnot be applied to the covariance matrix. Please specify the correction method (correction).")

    if(pivot == TRUE) {
      pivotvec            <- attr(covfactor, "pivot")
      covfactor           <- covfactor[, pivotvec]
    }

  } else if(correction == "qr") {

    if(pivot == TRUE) warning("This correction does not use pivoting method.")

    eigendecomp           <- eigen(covlocs)
    covfactor             <- tryCatch(qr.R(qr( diag(sqrt(eigendecomp$values)) %*% t(eigendecomp$vectors), tol = tol )), error = function(e) "Oops")

    if(identical(covfactor, "Oops")) stop("QR decomposition cannnot be applied. Please use another correction method (correction).")

  } else if(correction == "diag") {

    covfactor             <- tryCatch(chol(covlocs + diag(tol, n), pivot = pivot, tol = tol), error = function(e) "Oops")
    if(identical(covfactor, "Oops")) stop("Cholesky decomposition cannnot be applied to the modified covariance matrix. Please use larger tolerance value (tol).")

    if(pivot == TRUE) {
      pivotvec            <- attr(covfactor, "pivot")
      covfactor           <- covfactor[, pivotvec]
    }

  } else if(correction == "type-I") {

    diag(covlocs)         <- pmax(abs(diag(covlocs)), tol)
    covfactor             <- tryCatch(chol(covlocs, pivot = pivot, tol = tol), error = function(e) "Oops")
    if(identical(covfactor, "Oops")) stop("Cholesky decomposition cannnot be applied to the modified covariance matrix. Please use larger tolerance value (tol).")

    if(pivot == TRUE) {
      pivotvec            <- attr(covfactor, "pivot")
      covfactor           <- covfactor[, pivotvec]
    }

  } else if(correction == "type-II") {

    diag(covlocs)         <- pmax(diag(covlocs), tol)
    covfactor             <- tryCatch(chol(covlocs, pivot = pivot, tol = tol), error = function(e) "Oops")
    if(identical(covfactor, "Oops")) stop("Cholesky decomposition cannnot be applied to the modified covariance matrix. Please use larger tolerance value (tol).")

    if(pivot == TRUE) {
      pivotvec            <- attr(covfactor, "pivot")
      covfactor           <- covfactor[, pivotvec]
    }

  } else if(correction == "eigen-I") {

    if(pivot == TRUE) warning("This correction does not use pivoting method.")

    eigendecomp           <- eigen(covlocs)
    eigendecomp$values    <- pmax(abs(eigendecomp$values), tol) # ifelse(abs(eigendecomp$values) > tol, abs(eigendecomp$values), tol)
    covfactor             <- tryCatch(qr.R(qr( diag(sqrt(eigendecomp$values)) %*% t(eigendecomp$vectors) )), error = function(e) "Oops")

    if(identical(covfactor, "Oops")) stop("QR decomposition cannnot be applied. Please use larger tolerance value (tol).")

  } else if(correction == "eigen-II") {

    if(pivot == TRUE) warning("This correction does not use pivoting method.")

    eigendecomp           <- eigen(covlocs)
    eigendecomp$values    <- pmax(eigendecomp$values, tol) # ifelse(eigendecomp$values > tol, eigendecomp$values, tol)
    covfactor             <- tryCatch(qr.R(qr( diag(sqrt(eigendecomp$values)) %*% t(eigendecomp$vectors) )), error = function(e) "Oops")

    if(identical(covfactor, "Oops")) stop("QR decomposition cannnot be applied. Please use larger tolerance value (tol).")

  } else if(correction == "GMW81") {

    out <- tryCatch(.algorithm_GMW81(A = covlocs, beta = NULL, pivot = pivot, tol = tol), error = function(e) "Oops")

    if(identical(out, "Oops")) stop("GMW81 algorithm cannnot be applied. Please use larger tolerance value (tol).")

    if(pivot == FALSE) {
      covfactor <- sqrt(out$Dvec) * t(out$Lmat) # U = t(out$Lmat %*% diag(sqrt(out$Dvec)))
    } else if(pivot == TRUE) {
      covfactor <- sqrt(out$Dvec) * t(out$Lmat)
      covfactor <- covfactor[, out$Pvec]
    } else {
      stop("The argument pivot must be logical: TRUE or FALSE.")
    }

  } else if(correction == "SE99") {

    if(pivot == FALSE) warning("SE99 algorithm need to use pivoting method.")

    out <- tryCatch(.algorithm_SE99(A = covlocs, tol = tol), error = function(e) "Oops")

    if(identical(out, "Oops")) stop("SE99 algorithm cannnot be applied.")

    covfactor <- t(out$Lmat)
    covfactor <- covfactor[, out$Pvec]

  } else {

    stop("Not yet.")

  }

  covlocs.modified <- t(covfactor) %*% covfactor

  err.decomp <- sqrt(sum( (covlocs - covlocs.modified)^2 ))
  message("The decomposition error of the covariance matrix is ", err.decomp, " in the Frobenius norm.")

  y <- meanlocs + as.numeric(t(covfactor) %*% rnorm(n))

  return(list(y = y, meanlocs = meanlocs, covlocs = covlocs, covlocs.modified = covlocs.modified, err.decomp = err.decomp, correction = correction, tol = tol))
}

.algorithm_modchol <- function(A, delta, pivot) {

  n <- nrow(A)

  if(pivot == FALSE) {

    P <- NULL
    for(k in 1:(n-1)) {
      A[k, k] <- sqrt(A[k, k] + delta[k])
      A[(k+1):n, k] <- A[(k+1):n, k] / A[k, k]
      A[(k+1):n, (k+1):n] <- A[(k+1):n, (k+1):n] - tcrossprod(A[(k+1):n, k])
    }

    A[n, n] <- sqrt(A[n, n] + delta[n])
    A <- A * !upper.tri(A)

  } else {

    P <- seq(n)
    for(k in 1:(n-1)) {
      pvt             <- which.max(diag(A)[k:n])
      P[k:n]          <- c(P[k:n][pvt], P[k:n][-pvt])
      A               <- A[P, P]

      A[k, k] <- sqrt(A[k, k] + delta[k])
      A[(k+1):n, k] <- A[(k+1):n, k] / A[k, k]
      A[(k+1):n, (k+1):n] <- A[(k+1):n, (k+1):n] - tcrossprod(A[(k+1):n, k])
    }

    A[n, n] <- sqrt(A[n, n] + delta[n])
    A <- A * !upper.tri(A)
  }

  return(list(Lmat = A, Pvec = P))
}

.algorithm_modldl <- function(A, delta, pivot) {

  n <- nrow(A)

  L <- diag(n)
  D <- rep(NA, n)

  if(pivot == FALSE) {

    P     <- NULL
    A.k   <- A
    for(k in 1:(n-1)) {
      D[k]            <- A.k[1, 1] + delta[k]
      L[(k+1):n, k]   <- A.k[2:nrow(A.k), 1] / D[k]
      A.k             <- A.k[2:nrow(A.k), 2:nrow(A.k)] - tcrossprod(A.k[2:nrow(A.k), 1]) / D[k]
    }

    D[n] <- A.k[nrow(A.k), nrow(A.k)] + delta[n]

  } else {

    P     <- seq(n)
    A.k   <- A
    for(k in 1:(n-1)) {
      pvt             <- which.max(diag(A.k))
      P[k:n]          <- c(P[k:n][pvt], P[k:n][-pvt])
      A.k             <- A.k[c(pvt, (1:nrow(A.k))[-pvt]), c(pvt, (1:nrow(A.k))[-pvt])]
      L               <- L[P, P]

      D[k]            <- A.k[1, 1] + delta[k]
      L[(k+1):n, k]   <- A.k[2:nrow(A.k), 1] / D[k]
      A.k             <- A.k[2:nrow(A.k), 2:nrow(A.k)] - tcrossprod(A.k[2:nrow(A.k), 1]) / D[k]
    }

    D[n] <- A.k[nrow(A.k), nrow(A.k)] + delta[n]

  }

  return(list(Lmat = L, Dvec = D, Pvec = P))
}

.algorithm_GMW81 <- function(A, beta = NULL, pivot, tol = .Machine$double.eps) {

  n <- nrow(A)

  if(is.null(beta)) {
    eta <- max(abs( diag(A) ))
    xi <- max(abs( A * (1-diag(n)) ))
    beta <- sqrt(max( eta, xi / sqrt(n^2 - 1), .Machine$double.eps))
  }

  L   <- diag(n)
  D   <- rep(NA, n)
  P   <- seq(n)

  if(pivot == FALSE) {

    P     <- NULL
    A.k   <- A
    for(k in 1:(n-1)) {
      D[k]            <- max(tol, abs(A.k[1, 1]), max(A.k[2:nrow(A.k), 1]^2 / beta^2))
      L[(k+1):n, k]   <- A.k[2:nrow(A.k), 1] / D[k]
      A.k             <- A.k[2:nrow(A.k), 2:nrow(A.k)] - tcrossprod(A.k[2:nrow(A.k), 1]) / D[k]
    }

    D[n]  <- max(tol, abs(A.k[nrow(A.k), nrow(A.k)]))

  } else {

    P     <- seq(n)
    A.k   <- A
    for(k in 1:(n-1)) {
      pvt             <- which.max(diag(A.k))
      P[k:n]          <- c(P[k:n][pvt], P[k:n][-pvt])
      A.k             <- A.k[c(pvt, (1:nrow(A.k))[-pvt]), c(pvt, (1:nrow(A.k))[-pvt])]
      L               <- L[P, P]

      D[k]            <- max(tol, abs(A.k[1, 1]), max(A.k[2:nrow(A.k), 1]^2 / beta^2))
      L[(k+1):n, k]   <- A.k[2:nrow(A.k), 1] / D[k]
      A.k             <- A.k[2:nrow(A.k), 2:nrow(A.k)] - tcrossprod(A.k[2:nrow(A.k), 1]) / D[k]

    }

    D[n]  <- max(tol, abs(A.k[nrow(A.k), nrow(A.k)]))

  }

  return(list(beta = beta, Lmat = L, Dvec = D, Pvec = P))
}

.algorithm_bruteforceSE90 <- function(A, pivot) {

  n <- nrow(A)

  L <- diag(n)
  D <- pmax(diag(A), rowSums(abs(A)) - diag(abs(A)))

  if(pivot == FALSE) {

    P <- NULL
    A.k <- A
    for(k in 1:(n-1)) {
      L[(k+1):n, k]   <- A.k[2:nrow(A.k), 1] / D[k]
      A.k             <- A.k[2:nrow(A.k), 2:nrow(A.k)] - tcrossprod(A.k[2:nrow(A.k), 1]) / D[k]
    }

  } else {

    P     <- seq(n)
    A.k   <- A
    for(k in 1:(n-1)) {
      pvt             <- which.max(diag(A.k))
      P[k:n]          <- c(P[k:n][pvt], P[k:n][-pvt])
      A.k             <- A.k[c(pvt, (1:nrow(A.k))[-pvt]), c(pvt, (1:nrow(A.k))[-pvt])]
      L               <- L[P, P]
      D               <- D[P]

      L[(k+1):n, k]   <- A.k[2:nrow(A.k), 1] / D[k]
      A.k             <- A.k[2:nrow(A.k), 2:nrow(A.k)] - tcrossprod(A.k[2:nrow(A.k), 1]) / D[k]
    }

  }

  return(list(Lmat = L, Dvec = D, Pvec = P))
}

.algorithm_naiveSE90 <- function(A, pivot) {

  n <- nrow(A)

  L <- diag(n)
  D <- rep(NA, n)

  if(pivot == FALSE) {

    P <- NULL
    A.k <- A
    for(k in 1:(n-1)) {
      D[k]            <- max( A.k[1, 1], sum(abs(A.k[2:nrow(A.k), 1])) )
      L[(k+1):n, k]   <- A.k[2:nrow(A.k), 1] / D[k]
      A.k             <- A.k[2:nrow(A.k), 2:nrow(A.k)] - tcrossprod(A.k[2:nrow(A.k), 1]) / D[k]
    }

    D[n] <- A.k[nrow(A.k), nrow(A.k)]

  } else {

    P     <- seq(n)
    A.k   <- A
    for(k in 1:(n-1)) {
      pvt             <- which.max(diag(A.k))
      P[k:n]          <- c(P[k:n][pvt], P[k:n][-pvt])
      A.k             <- A.k[c(pvt, (1:nrow(A.k))[-pvt]), c(pvt, (1:nrow(A.k))[-pvt])]
      L               <- L[P, P]

      D[k]            <- max( A.k[1, 1], sum(abs(A.k[2:nrow(A.k), 1])) )
      L[(k+1):n, k]   <- A.k[2:nrow(A.k), 1] / D[k]
      A.k             <- A.k[2:nrow(A.k), 2:nrow(A.k)] - tcrossprod(A.k[2:nrow(A.k), 1]) / D[k]
    }

    D[n]  <- A.k[nrow(A.k), nrow(A.k)]

  }

  return(list(Lmat = L, Dvec = D, Pvec = P))
}

.algorithm_SE99 <- function(A, tol) { # always use pivoting strategy

  outmodchol <- mvLSW:::modchol(A, tol = tol)

  return(list(Lmat = outmodchol$L, Pvec = outmodchol$P))

}

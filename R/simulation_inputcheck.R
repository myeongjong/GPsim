####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions for
###
###   Description:
###
###
####################################################################################


.input_check_locs_n_p <- function(locs, n, p) {

  if( !is.null(n) ) {
    if(class(n) != 'numeric') stop("The number of locations (n) is not numeric.")
    if(round(n) != n) stop("The number of locations (n) is not an integer.")
    if(n <= 0) stop("The number of locations (n) is not a positive integer.")
    if(n == 1) stop("The number of locations (n) must be larger than 1.")
  }

  if( !is.null(p) ) {
    if(class(p) != 'numeric') stop("The dimension of domain (p) is not numeric.")
    if(round(p) != p) stop("The dimension of domain (p) is not an integer.")
    if(p <= 0) stop("The dimension of domain (p) is not a positive integer.")
  }

  if(is.null(locs)) {

    if(is.null(n) | is.null(p)) {
      stop("You must either provide locations (locs) or specify the number of locations (n), the dimension of domain (p), and how to select locations (locs).")
    } else {
      message("Warning: You did not specify how to select locations (locs). Locations are randomly selected.")
      locs    <- as.data.frame(matrix(runif(n * p, 0, 1), n, p))
    }

  } else if(is.matrix(locs)) {

    if(is.null(n) & is.null(p)) {
      message("Note: Both the number of locations (n) and the dimension of domain (p) are defined with respect to the given locations (locs).")
      n <- nrow(locs)
      p <- ncol(locs)
    } else if(is.null(n) & !is.null(p)) {
      if(ncol(locs) == p) {
        message("Note: The number of locations (n) is defined with respect to the given locations (locs).")
        n <- nrow(locs)
      } else {
        stop("The dimension of domain (p) is not compatible with the locations (locs).")
      }
    } else if(!is.null(n) & is.null(p)) {
      if(nrow(locs) == n) {
        message("Note: The dimension of domain (p) is defined with respect to the given locations (locs).")
        p <- ncol(locs)
      } else {
        stop("The number of locations (n) is not compatible with the locations (locs).")
      }
    } else if(!( is.null(n) | is.null(p))){
      if(identical(as.numeric(dim(locs)), c(n, p))) {
        message("Note: Neither the number of locations (n) nor the dimension of domain (p) are necessary.")
      } else {
        stop("The locations (locs) are not compatible with either the number of locations (n) or the dimension of domain (p).")
      }
    }

    locs <- as.data.frame(locs)

  } else if(is.data.frame(locs)) {

    if(is.null(n) & is.null(p)) {
      message("Note: Both the number of locations (n) and the dimension of domain (p) are defined with respect to the given locations (locs).")
      n <- nrow(locs)
      p <- ncol(locs)
    } else if(is.null(n) & !is.null(p)) {
      if(ncol(locs) == p) {
        message("Note: The number of locations (n) is defined with respect to the given locations (locs).")
        n <- nrow(locs)
      } else {
        stop("The dimension of domain (p) is not compatible with the locations (locs).")
      }
    } else if(!is.null(n) & is.null(p)) {
      if(nrow(locs) == n) {
        message("Note: The dimension of domain (p) is defined with respect to the given locations (locs).")
        p <- ncol(locs)
      } else {
        stop("The number of locations (n) is not compatible with the locations (locs).")
      }
    } else if(!( is.null(n) | is.null(p))){
      if(identical(as.numeric(dim(locs)), c(n, p))) {
        message("Note: Neither the number of locations (n) nor the dimension of domain (p) are necessary.")
      } else {
        stop("The locations (locs) are not compatible with either the number of locations (n) or the dimension of domain (p).")
      }
    }

  } else if(locs == "grid") {

    if(is.null(n) | is.null(p)) {
      stop("You must provide the number of locations (n) and the dimension of domain (p) when specifying how to select locations (locs).")
    } else {
      grid.oneside    <- seq(from = 0, to = 1, length.out = n^(1/p))
      locs            <- expand.grid( as.data.frame(matrix(grid.oneside, nrow = length(grid.oneside), ncol = p, byrow = F)) )

      if(n != nrow(locs)) {
        n             <- nrow(locs)
        message("Warning: The arguments n and p are not compatible. n will be ", n, ".")
      }
    }

  } else if(locs == "random") {

    if(is.null(n) | is.null(p)) {
      stop("You must provide the number of locations (n) and the dimension of domain (p) when specifying how to select locations (locs).")
    } else {
      locs            <- as.data.frame(matrix(runif(n * p, 0, 1), n, p))
    }

  } else {
    stop("The input locs is not valid. Please refer the description.")
  }

  colnames(locs) <- paste0('s', seq(p))
  return(list(locs = locs, n = n, p = p))
}

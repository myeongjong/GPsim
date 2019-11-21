####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions for simulating Gaussian Processes (GPs).
###
###   Description:
###
###
####################################################################################

#' Title
#'
#' @param temp1
#' @param temp2
#'
#' @return
#' @export
#'
#' @examples
simulate_gp <- function(temp1, temp2) {
  return(temp1 + temp2)
}

.simulate_gp_stationary <- function() {
  return(1)
}

.simulate_gp_nonstationary <- function() {
  return(1)
}

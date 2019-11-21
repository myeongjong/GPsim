####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions for implementing baseline models.
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
corrvecchia <- function(temp1, temp2) {
  return(temp1 + temp2)
}

.corrvecchia_knownCovparams <- function() {
  return(1)
}

####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions for
###
###   Description:
###
###
####################################################################################


###  function to compute KL divergence
kldiv <- function(covmat0,covmat1,mu0=rep(0,nrow(covmat0)),mu1=rep(0,nrow(covmat0)),chol1=NULL){

  # multivariate normal KL divergence between
  # true distribution N(mu0,covmat0) and approx distr N(mu1,covmat1)
  #   if chol of covmat1 has already been computed, provide chol1 instead of covmat1

  n <- nrow(covmat0)
  chol0 <- t(chol(covmat0))
  if(is.null(chol1)) chol1 <- t(chol(covmat1))

  # trace term
  M <-   solve( chol1, chol0 )
  traceterm <- sum( M^2 ) - n

  # log det term
  logdetterm <- 2*sum(log(diag(chol1))) - 2*sum(log(diag(chol0)))

  # mean term (equal to zero if means are zero)
  meandiff=mu1-mu0
  # t(meandiff)%*%solve(covmat1,meandiff)
  temp=forwardsolve(chol1,meandiff)
  meanterm=sum(temp^2)

  # kldiv
  kldiv <- 1/2*( traceterm + logdetterm + meanterm )
  return(kldiv)

}

####################################################################################
###   Author: Myeongjong Kang (kmj.stat@gmail.com)
###
###   Overview: This script includes several R functions for
###
###   Description:
###
###
####################################################################################


.order_coordinate <- function (locs, coordinate) {
  if(ncol(locs) == 1){
    order(rowSums(locs[, drop = FALSE]))
  } else {
    order(rowSums(locs[, coordinate, drop = FALSE]))
  }
}


.order_maxmin_euclidean <- function(locs)
{
  n     <- nrow(locs)
  p     <- ncol(locs)
  ord   <- rep(NA, n)
  cen   <- t(as.matrix(colMeans(locs)))

  ord[1]        <- which.min(rowSums((locs - matrix(as.numeric(cen), nrow = n, ncol = p, byrow = T))^2))
  cand.argmax   <- seq(n)[seq(n) != ord[1]]

  cdist         <- fields::rdist(locs[cand.argmax, ], matrix(locs[ord[1], ], nrow = 1, ncol = 2))
  ord[2]        <- cand.argmax[which.max(as.numeric(cdist))]
  cand.argmax   <- cand.argmax[cand.argmax != ord[2]]

  for(j in 3:(n-1)){
    cdist         <- fields::rdist(locs[cand.argmax, ], locs[ord[seq(j-1)], ])
    cdist         <- Rfast::rowMins(cdist, value = T)
    ord[j]        <- cand.argmax[which.max(cdist)]
    cand.argmax   <- cand.argmax[cand.argmax != ord[j]]
  }

  ord[n]        <- cand.argmax

  return(ord)
}


### CAUTION: This function can cause numerical issue. Please use the 'correlation()' function, instead. ###
.distance_correlation <- function(locs, covmodel, covparms, def.dist)
{
  if(is.function(covmodel)) {
    covparms[1]   <- 1 # correlation function
    dist.matrix   <- if(is.null(def.dist)) 1 - covmodel(locs, covparms) else 1 - abs(covmodel(locs, covparms)) # 1-rho
  } else if(is.matrix(covmodel)) {
    dist.matrix   <- if(is.null(def.dist)) 1 - covmodel / covparms[1] else 1 - abs(covmodel) / covparms[1]
  } else {
    dist.matrix   <- 1 - diag(1, nrow = nrow(locs), ncol = nrow(locs))
  }

  return(dist.matrix)
}

.correlation <- function(locs, covmodel, covparms, def.dist)
{
  if(is.function(covmodel)) {
    covparms[1]   <- 1 # correlation function
    corr.matrix   <- if(is.null(def.dist)) covmodel(locs, covparms) else abs(covmodel(locs, covparms)) # 1-rho
  } else if(is.matrix(covmodel)) {
    corr.matrix   <- if(is.null(def.dist)) covmodel / covparms[1] else abs(covmodel) / covparms[1]
  } else {
    corr.matrix   <- diag(1, nrow = nrow(locs), ncol = nrow(locs))
  }

  return(corr.matrix)
}


.order_maxmin_correlation <- function(locs, dinv, covmodel, covparms, initial.pt = NULL)
{
  n     <- nrow(locs)
  p     <- ncol(locs)
  ord   <- rep(NA, n)

  if( is.null(initial.pt) ){
    ord[1]        <-  which.max(rowSums(dinv))
  } else if( is.numeric(initial.pt) ) {
    ord[1]        <- initial.pt
  } else if( initial.pt == 'center' ) {
    cen           <- t(as.matrix(colMeans(locs)))
    ord[1]        <- which.min(rowSums((locs - matrix(as.numeric(cen), nrow = n, ncol = p, byrow = T))^2))
  } else if( initial.pt == 'random') { # at random
    ord[1]        <- sample(1:n, 1)
  } else {
    ord[1]        <-  which.max(rowSums(dinv))
  }
  cand.argmax   <- seq(n)[seq(n) != ord[1]]

  ind           <- as.matrix(expand.grid(cand.argmax, ord[1]))
  cdist         <- dinv[ind]
  ord[2]        <- cand.argmax[which.min(as.numeric(cdist))]
  cand.argmax   <- cand.argmax[cand.argmax != ord[2]]

  for(j in 3:(n-1)){
    ind           <- as.matrix(expand.grid(cand.argmax, ord[seq(j-1)]))
    cdist         <- matrix(dinv[ind], nrow = length(cand.argmax), ncol = j-1, byrow = F)
    cdist         <- Rfast::rowMaxs(cdist, value = T)
    ord[j]        <- cand.argmax[which.min(cdist)]
    cand.argmax   <- cand.argmax[cand.argmax != ord[j]]
  }

  ord[n]        <- cand.argmax

  return(ord)
}


.conditioning_nn <- function(m, dist.matrix)
{
  # initialize an output object NN which is a n*n matrix
  n     <- nrow(dist.matrix)
  NN    <- matrix(rep(NA, n * (m + 1)), nrow = n, ncol = m + 1) ; NN[1, 1] <- 1

  # Find the nearest neighbor conditioning set for each i-th location using the 'dist_to_knn()' function
  for(i in 2:n) {
    k               <- min(i, m + 1) # the number of neighbors of the i-th observation
    NN[i, seq(k)]   <- scanstatistics::dist_to_knn(dist.matrix[seq(i), seq(i)], k = k)[i, seq(k)]
  }

  return(NN)
}


.U_sparsity <- function(locs, NNarray, obs, Cond) {
  # locs are the locations
  # row i of NNarray gives the neighbors of location i
  # i-th element of observed is TRUE if loc i is observed

  nnp   <- nrow(locs) # total number of locs (incl unobserved)
  n     <- sum(obs) # number of obs locs
  size  <- nnp + n # size (i.e., number of rows and columns in U)

  # Calculate how many nonzero entries we need to compute
  nentries <- sum( !is.na(NNarray) )

  # latent_map specifies which row in U corresponds to which latent variable
  # observed_map specifies which row in U corresponds to which obs variable
  cur <- 1
  latent_map <- rep(NA, n)
  observed_map <- rep(NA, n)
  for(k in 1:nnp){
    latent_map[k] <- cur
    cur <- cur+1
    if( obs[k] ){
      observed_map[k] <- cur
      cur <- cur+1
    }
  }

  # pre-reorder everything
  revNNarray <- t(apply(NNarray, 1, rev)) # column-reversed NNarray
  revCondOnLatent <- t(apply(Cond, 1, rev)) # column-reversed CondOnLatent T/F

  # calculate the indexing of latent variables in U
  rowpointers = colindices = rep(NA, nentries)
  rowpointers[1] = colindices[1] = 1
  cur <- 0
  for( k in 1:nnp){
    inds    <- revNNarray[k,]
    inds0   <- inds[!is.na(inds)]
    n0 <- length(inds0)
    revCond <- revCondOnLatent[k,!is.na(inds)]

    cur_row <- latent_map[k]
    cur_cols <- rep(NA, n0)
    cur_cols[revCond]  <- latent_map[inds0[revCond]]
    cur_cols[!revCond] <- observed_map[inds0[!revCond]]

    # store pointers to appropriate rows and columns
    if( k > 1 ){
      rowpointers[cur + (1:n0)] <- as.integer(rep(cur_row, n0))
      colindices[cur + (1:n0)] <- as.integer(cur_cols)
    }
    cur <- cur+n0
  }

  # separately, calculate indexing of obs variables in U
  Zrowpointers = Zcolindices = rep(NA, 2 * n)
  cur <- 0
  for(k in 1:nnp){
    if( obs[k] ){
      cur_row <- observed_map[k]
      cur_cols <- c( latent_map[k], observed_map[k] )
      Zrowpointers[cur + (1:2)] <- as.integer(rep(cur_row, 2))
      Zcolindices[cur + (1:2)] <- as.integer(cur_cols)
      cur <- cur + 2
    }
  }

  # combine pointers for latent and obs variables
  allrowpointers  <- c(rowpointers, Zrowpointers)
  allcolindices   <- c(colindices, Zcolindices)

  # number of cores to be used for creating U later
  n.cores <- parallel::detectCores(all.tests = FALSE, logical = TRUE)

  return(list(revNNarray = revNNarray, revCond = revCondOnLatent, n.cores = n.cores, size = size, rowpointers = allrowpointers, colindices = allcolindices, y.ind = latent_map, observed_map = observed_map))

}

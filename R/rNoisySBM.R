#' Simulate Noisy SBM
#'
#' \code{rNoisyMBM} simulates a collection of networks which are noisy versions of un underlying network.  This network may be directed or not. See vignette for more informations
#' @param  nbNodes        : number of nodes
#' @param  directed      : directed or not (directed = TRUE or FALSE)
#' @param  blockProp     : proportions of the Blocks. Vector of size nbBlocks
#' @param  connectParam  : probabilities of connections inside and between blocks. Matrix of size nbBlocks.
#' @param  emissionParam : parameters of the emission of the Noisy SBM. List of two terms : noEdgeParam and edgeParam. Each element is a list of means (vector or size nbScores) and variance matrix (square matrix of size nbScores)
#' @param  seed          : set the seed for the random simulation (default value  = NULL)
#' @return A list of containing the noisy networks (noisyNetworks a list of length nbScores of matrices of dimension nbNodes x nbNodes), the true underlying network (trueNetwork) and the clustering of the nodes (memberships)
#' @examples
#' nbNodes  <- 100
#' directed <- TRUE
#' blockProp <- c(1/3,1/2,1/6)
#' nbBlocks   <- length(blockProp)
#' connectParam <- matrix(rbeta(nbBlocks^2,1.5,1.5 ),nbBlocks,nbBlocks)
#' connectParam <- 0.5*(connectParam + t(connectParam))
#' emissionParam <- list()
#' nbScores <- 4
#' emissionParam$noEdgeParam <- list(mean=rep(0,nbScores));
#' emissionParam$noEdgeParam$var <- diag(0.1,nrow = nbScores,ncol = nbScores)
#' emissionParam$EdgeParam <- list( mean= 1:nbScores)
#' emissionParam$EdgeParam$var <-  diag(0.1,nrow = nbScores,ncol = nbScores)
#' dataSim <- rNoisySBM(nbNodes,directed <- TRUE, blockProp,connectParam,emissionParam,seed = NULL)
#'
#' @export
rNoisySBM = function(nbNodes,
  directed,
  blockProp,
  connectParam,
  emissionParam,
  seed = NULL) {
  #---------------------  TESTS on dimensions

  #--- on the SBM parameters
  K <- length(blockProp)
  Krow <- nrow(connectParam)
  Kcol <- ncol(connectParam)
  if ((K != Krow) |
      (K !=  Kcol)) {
    stop('Non matching dimensions between blockProp and connectParam.')
  }
  if ( directed  &
      (!isSymmetric(connectParam))) {
    stop('For a non-directed (symetric) network connectParam should be a symetric matrix ')
  }

  #-- on the means of the emission parameters

  dE <- length(emissionParam$EdgeParam$mean)
  dnoE <- length(emissionParam$noEdgeParam$mean)
  if (dE != dnoE) {
    stop('The two means of the Emission distributions should be of the same sizes')
  }
  d <- dE


  if (d > 1) {
    if (!isSymmetric(emissionParam$EdgeParam$var) |
        !isSymmetric(emissionParam$noEdgeParam$var)) {
      stop('One (or both) of the var parameters of the emission distribution are not Symetric')
    }
    if (d != nrow(emissionParam$EdgeParam$var) |
        d != nrow(emissionParam$noEdgeParam$var)) {
      stop(' Check the sizes od the eimssion variance matrices')
    }
  }




  #---- SEED------------------------
  set.seed(seed)

  Z <- stats::rmultinom(nbNodes, 1, blockProp)

  Pconnect <- t(Z) %*% connectParam  %*% Z
  G <- matrix(stats::rbinom(nbNodes ^ 2, 1, Pconnect), nbNodes, nbNodes)
  GArray <- replicate(d, G, simplify = "array")
  X1 <-
    mvtnorm::rmvnorm(nbNodes ^ 2,
      emissionParam$EdgeParam$mean,
      emissionParam$EdgeParam$var)
  X0 <-
    mvtnorm::rmvnorm(nbNodes ^ 2,
      emissionParam$noEdgeParam$mean,
      emissionParam$noEdgeParam$var)
  X1Array <- array(X1, c(nbNodes, nbNodes, d))
  X0Array <- array(X0, c(nbNodes, nbNodes, d))
  XArray <- X1Array * GArray + X0Array * (1 - GArray)
  X <- lapply(1:d, function(i) {U <- XArray[, , i]; diag(U) <- NA; return(U)})
  if (directed){
    where <- upper.tri(G)
    G[where] <- t(G)[where]
    diag(G) <- rep(0,nbNodes)
    X <- lapply(1:d, function(i) {
      0.5 * (X[[i]] + t(X[[i]]))
    })
  }

  res = list(
    trueNetwork = G,
    noisyNetworks = X,
    memberships = t(Z)
  )
  return(res)

}

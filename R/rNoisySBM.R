#' Simulate Noisy SBM
#'
#' \code{rNoisyMBM} simulates a collection of networks which are noisy versions of un underlying network.  This network may be directed or not. See vignette for more informations
#' @param  nNodes        : number of nodes
#' @param  symmetric     : symmetric (directed or not) (symmetric = TRUE or FALSE)
#' @param  mixtureParam  : proportions of the Blocks. Vector of size K
#' @param  connectParam  : probabilities of connections inside and between blocks. Matrix of size K.
#' @param  emissionParam : parameters of the emission of the Noisy SBM. List of two terms : noEdgeParam and edgeParam. Each element is a list of means (vector or size d) and variance matrix (square matrix of size d)
#' @param  seed          : set the seed for the random simulation (default value  = NULL)
#' @return A list of lists containing the noisy networks (noisyNetworks), the true underlying network trueNetwork and
#'         Each element of  list_net corresponds to a network : each network is a list containing  the matrix (mat) , the type of network(diradj, adj, inc), the functional group in row (rowFG) and the functional group in columns (colFG)
#' @examples
#' nNodes  = 100
#' directed = TRUE
#' mixtureParam = c(1/3,1/2,1/6)
#' nBlocks   = length(mixtureParam)
#' connectParam <- matrix(rbeta(nBlocks^2,1.5,1.5 ),nBlocks,nBlocks)
#' connectParam <- 0.5*(connectParam + t(connectParam))
#' emissionParam <- list()
#' d <- 4
#' emissionParam$noEdgeParam = list(mean=rep(0,d),var = diag(0.1,nrow = d,ncol = d))
#' emissionParam$EdgeParam = list( mean= 1:d,var =  diag(0.1,nrow = d,ncol = d))
#' dataSim <- rNoisySBM(nNodes,symmetric = TRUE, mixtureParam,connectParam,emissionParam,seed = NULL)
#'
#' @export
rNoisySBM = function(nNodes,
  symmetric,
  mixtureParam,
  connectParam,
  emissionParam,
  seed = NULL) {
  #---------------------  TESTS on dimensions

  #--- on the SBM parameters
  K <- length(mixtureParam)
  Krow <- nrow(connectParam)
  Kcol <- ncol(connectParam)
  if ((K != Krow) |
      (K !=  Kcol)) {
    stop('Non matching dimensions between mixtureParam and connectParam.')
  }
  if ( symmetric  &
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

  Z <- stats::rmultinom(nNodes, 1, mixtureParam)

  Pconnect <- t(Z) %*% connectParam  %*% Z
  G <- matrix(stats::rbinom(nNodes ^ 2, 1, Pconnect), nNodes, nNodes)
  GArray <- replicate(d, G, simplify = "array")
  X1 <-
    mvtnorm::rmvnorm(nNodes ^ 2,
      emissionParam$EdgeParam$mean,
      emissionParam$EdgeParam$var)
  X0 <-
    mvtnorm::rmvnorm(nNodes ^ 2,
      emissionParam$noEdgeParam$mean,
      emissionParam$noEdgeParam$var)
  X1Array <- array(X1, c(nNodes, nNodes, d))
  X0Array <- array(X0, c(nNodes, nNodes, d))
  XArray <- X1Array * GArray + X0Array * (1 - GArray)
  X <- lapply(1:d, function(i) {
    XArray[, , i]
  })
  if (symmetric){
    where <- upper.tri(G)
    G[where] <- t(G)[where]
    diag(G) <- rep(0,nNodes)
    X <- lapply(1:d, function(i) {
      0.5 * (X[[i]] + t(X[[i]]))
    })
  }

  res = list(
    trueNetwork = G,
    noisyNetworks = X,
    memberships = Z
  )
  return(res)

}

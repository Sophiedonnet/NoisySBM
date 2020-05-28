#-------------------------------------------------------------------------------
# Computes the number of (symmetric) pairs
#-------------------------------------------------------------------------------
#' n2nbPairs : compute the number of Dyads from the number of nodes for a network
#' @param n : number of nodes
#' @param symmetric  : TRUE is the network is not  directed. FALSE otherwise
#' @param diag  : FALSE if the self interactions are neglected. TRUE otherwise (default value =  FALSE)
#' @return number of dyades
#' @export
#' @examples
#' n <- 10;
#' N <- n2nbPairs(n, symmetric = TRUE, diag=FALSE)
n2nbPairs <- function(n, symmetric, diag = FALSE){

  if (symmetric) {
    N <- ifelse(diag, n*(n + 1)/2, n*(n - 1)/2)
  }else{
    N <- ifelse(diag, N <- n^2, n*(n - 1))
  }
  return(N)
}

#-------------------------------------------------------------------------------
# Computes the individual from the number of (symmetric) pairs
#-------------------------------------------------------------------------------
#' nbPairs2n : compute the number of nodes of a network from its number of dyads
#' @param N  : number of dyads
#' @param symmetric  : TRUE is the network is not  directed. FALSE otherwise
#' @param diag  : FALSE if the self interactions are neglected. TRUE otherwise (default value =  FALSE)
#' @return number of dyades
#' @export
#' @examples
#' N <- 45;
#' n <- nbPairs2n(n, symmetric = TRUE, diag=FALSE)
nbPairs2n <- function(N, symmetric, diag=FALSE){
  if (symmetric) {
    n <- ifelse(diag, (-1 + sqrt(1 + 8*N))/2 , (1 + sqrt(1 + 8 * N))/2)
  }else{
    n <- ifelse(diag, sqrt(N) , (1 + sqrt(1 + 4*N))/2)
  }
  return(n)
}

#-------------------------------------------------------------------------------
# Builds a (symmetric) matrix from a vector
#-------------------------------------------------------------------------------
#' vect2Mat : transform a vector into a  matrix  (replicating symetric elements if needed)
#'
#' @param V : vector
#' @param symmetric  : TRUE is the network is not  directed. FALSE otherwise
#' @param diag  : FALSE if the self interactions are neglected. TRUE otherwise (default value =  FALSE)
#' @return matrix fulfilled with the elements of the matrix
#' @export
#' @examples
#' V <- c(0,15);
#' M <- vect2Mat(V, symmetric = TRUE, diag=FALSE)
vect2Mat <- function(V, symmetric, diag=FALSE)
{

  N <- length(V); n <- nbPairs2n(N, symmetric=symmetric, diag=diag);
  M <- matrix(0,n,n)

  if (symmetric == TRUE){
    # n <- ifelse(diag, (-1 + sqrt(1 + 8*N))/2 , (1 + sqrt(1 + 8 * N))/2)
    # M <- matrix(0,n,n)
    M[lower.tri(M, diag = diag)] <- V
    B <- t(M)
    diag(B) <- 0
    M <- B + M
  }else{
    # n <- ifelse(diag, sqrt(N) , (1 + sqrt(1 + 4*N))/2)
    # M <- matrix(0,n,n)
    where <- (lower.tri(M, diag = diag) +  upper.tri(M, diag = diag)) > 0
    M[where] <- V
  }
  return(M)
}

#-------------------------------------------------------------------------------
# Builds a vector from a (symmetric) matrix
#-------------------------------------------------------------------------------
#' mat2Vect : Builds a vector from a (symmetric) matrix
#'
#' @param V : vector
#' @param symmetric  : TRUE is the network is not  directed. FALSE otherwise
#' @param diag  : FALSE if the self interactions are neglected. TRUE otherwise (default value =  FALSE)
#' @return matrix fulfilled with the elements of the matrix
#' @export
#' @examples
#' V <- c(0,15);
#' M <- vect2Mat(V, symmetric = TRUE, diag=FALSE)
#' V2 <- mat2Vect(M, symmetrice = TRUE, diag = FALSE)
mat2Vect <- function(M, symmetric, diag=FALSE)
{
  if (symmetric) {
    where <- lower.tri(M, diag = diag)
  } else {
    where <- (lower.tri(M, diag=diag) +  upper.tri(M, diag = diag)) > 0
  }
  V <- M[where]
  return(V)
}

#-------------------------------------------------------------------------------
#- Get the indices corresponding to a list of (symmetric) pairs
#-------------------------------------------------------------------------------
indices <- function(n, symmetric, diag=FALSE)
{
  N <- n2nbPairs(n, symmetric = symmetric, diag=diag)
  S <- vect2Mat(1:N, symmetric = symmetric, diag=diag)
  if (symmetric) {S[upper.tri(S)] = 0}
  res <- which(S!=0, arr.ind=TRUE)
  return(res)
}



#------------------------------------------------------------------------------
#---- Transorm the list of the matrices into a unique matrix
#------------------------------------------------------------------------------
#' scoreList2scoreMat : transform a list of Score matrices (matrices of size nxn) into a matrix
#' @param listScores : list of scores vector
#' @param symmetric  : TRUE is the network is not  directed. FALSE otherwise
#' @return one matrix with 'nbDyads' rows and d columns
#' @export
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
#' emissionParam$edgeParam <- list( mean= 1:nbScores)
#' emissionParam$edgeParam$var <-  diag(0.1,nrow = nbScores,ncol = nbScores)
#' dataSim <- rNoisySBM(nbNodes,directed = TRUE, blockProp,connectParam,emissionParam,seed = NULL)
#' S <- scoreList2scoreMat(dataSim$noisyNetworks , symmetric = FALSE)
scoreList2scoreMat <- function(listScores,symmetric){

  S <- listScores
  d <- length(S);
  #n <- nrow(S[[1]])
  mat_S <- sapply(1:d, function(q) {mat2Vect(S[[q]], symmetric = symmetric, diag = F)})
  # of dimension N  * d where N = n(n-1)/2 if symmetric, n(n-1) otherwise
  return(mat_S)
}





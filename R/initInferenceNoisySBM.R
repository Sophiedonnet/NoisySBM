#' Initialization of the inference procedure
#'
#' This function initialises the inference method by mixing a Gaussian mixture on the scores of each pair of nodes and a SBM on the resulting estimated network G
#' @param scoreList a list of the Scores (matrices of size nbNodes x nbNodes)
#' @param directed  a logical : TRUE if the underlying network is directed,  FALSE otherwise (default value FALSE).
#'
#' @return An initialisation point (list, see Details here after) for the VEM algortihms for several numbers of blocks
#' @details The output is a list containing the following entries.
#'  \itemize{
#'  \item{"G": }{estimated underlyling network (binary)}
#'  \item{"psi": }{estimated probability for each pair of nodes to be linked or not}
#'  \item{"tau": }{list of posterior probabilities for each not to be in a blocks. The length of the list corresponds to the various models explored by the SBM }
#'  \item{"ICL": }{Penalized likelihood for each of the estimated model}
#'  \item{"eta": }{Additional quantites used in the VEM}
#' }

#' @examples
#' nbNodes  = 60
#' directed = TRUE
#' blockProp = c(1/3,1/2,1/6)
#' nbBlocks   = length(blockProp)
#' connectParam <- matrix(rbeta(nbBlocks^2,1.5,1.5 ),nbBlocks,nbBlocks)
#' connectParam <- 0.5*(connectParam + t(connectParam))
#' emissionParam <- list()
#' nbScores <- 2
#' emissionParam$noEdgeParam <- list(mean=rep(0,nbScores));
#' emissionParam$noEdgeParam$var <- diag(0.1,nrow = nbScores,ncol = nbScores)
#' emissionParam$EdgeParam <- list( mean= 1:nbScores)
#' emissionParam$EdgeParam$var <-  diag(0.1,nrow = nbScores,ncol = nbScores)
#' dataSim <- rNoisySBM(nbNodes,directed, blockProp,connectParam,emissionParam,seed = NULL)
#' init <- initInferenceNoisySBM(dataSim$noisyNetworks, directed = FALSE)
#' @importFrom mclust mclustBIC Mclust
#' @export
#'
initInferenceNoisySBM <- function(scoreList, directed = FALSE){

  nbScores <- length(scoreList)
  scoreMat <- sapply(1:nbScores , function(q) {mat2Vect(scoreList[[q]], symmetric = !directed, diag = F)})

  N <- nrow(scoreMat)
  #--------  Mclust on the d scores dyad per dyad
  param_gm <- mclust::Mclust(scoreMat, G = 2, verbose = FALSE)
  psi <- param_gm$z
  G <- param_gm$classification - 1
  mu  <- param_gm$parameters$mean
  test_G <- rowMeans(t(mu)) #identify G = 0  and G =1
  if (test_G[1] > test_G[2]) {
    psi <- psi[,c(2,1)]
    G = 1 - G
  }


  #------------------ init of SBM parameters
  membership_type <- ifelse(directed, "SBM", "SBM_sym")
  param_sbm <- blockmodels::BM_bernoulli(membership_type, adj = vect2Mat(G, symmetric = !directed),
                                         plotting = '',
                                         verbosity = 0)
  param_sbm$estimate()
  Kmax <- length(param_sbm$memberships)
  tau_init <-  lapply(1:Kmax, function(K){param_sbm$memberships[[K]]$Z})
  eta  <-  lapply(1:Kmax,  function(K){array(rep(psi[, 1], K * K), c(N, K, K))})
  #-------------------------------------------------------------------------------

  res <- list(psi = psi, tau = tau_init, eta = eta, ICL = param_sbm$ICL, G = G)
  return(res)

}


#-------------------------------------------------------------------------------
# Initialization of the inference procedure
#-------------------------------------------------------------------------------
initInference <- function(mat_S, directed){

  N <- nrow(mat_S)
  #--------  Mclust on the d scores dyad per dyad
  param_gm <- mclust::Mclust(mat_S, G = 2, verbose = FALSE)
  Psi <- param_gm$z
  G <- param_gm$classification - 1
  mu  <- param_gm$parameters$mean
  test_G <- rowMeans(t(mu)) #identify G = 0  and G =1
  if (test_G[1] > test_G[2]) {
    Psi <- Psi[,c(2,1)]
    G = 1 - G
  }


  #------------------ init of SBM parameters
  membership_type <- ifelse(directed, "SBM", "SBM_sym")
  param_sbm <- blockmodels::BM_bernoulli(membership_type, adj = vect2Mat(G, symmetric=!directed),
                                       plotting = '',
                                       verbosity = 0)
  param_sbm$estimate()
  Kmax <- length(param_sbm$memberships)
  tau_init <-  lapply(1:Kmax, function(K){param_sbm$memberships[[K]]$Z})
  eta0 <-  lapply(1:Kmax,  function(K){array(rep(Psi[, 1], K * K), c(N, K, K))})
  #-------------------------------------------------------------------------------

  res <- list(Psi = Psi, tau = tau_init, eta0 = eta0, ICL = param_sbm$ICL, G = G)
  return(res)

}


#-------------------------------------------------------------------------------
# M step of the VEM algorithm
#-------------------------------------------------------------------------------
mStep <- function(scoreMat, psi, tau, eta, directed){

  # scoreMat <- mat_S; psi <- initQuanti$Psi; directed <- FALSE
  # nBlocks <- 3; tau <- initQuanti$tau[[nBlocks]]; eta <- initQuanti$eta0[[nBlocks]];

  # Dimensions
  N <- nrow(scoreMat); d <- ncol(scoreMat); nBlocks <- ncol(tau)

  # Proportions
  blockProp <- colMeans(tau)

  # Connection probabilities
  connectParam <- matrix(0, nBlocks, nBlocks)
  sapply(1:nBlocks, function(k){sapply(1:nBlocks, function(l){
    tauVec <- mat2Vect(tau[, k]%o%tau[, l], symmetric=!directed, diag=FALSE)
    connectParam[k, l] <<- tauVec%*%eta[, k, l] / sum(tauVec)
  })})

  # Emission distributions: mu and Sigma
  mu <- matrix(0, 2, d); Sigma <- array(dim=c(2, d, d))
  sapply(1:2, function(g){
    mu[g, ] <<- t(psi[, g]) %*% scoreMat / sum(psi[, g])
    Sigma[g, , ] <<- t(scoreMat) %*% diag(psi[, g]) %*% scoreMat / sum(psi[, g])
    Sigma[g, , ] <<- Sigma[g, , ] - mu[g, ]%o%mu[g, ]
    Sigma[g, , ] <<- .5*(Sigma[g, , ] + t(Sigma[g, , ]))
  })
  emissionParam <- list(noEdgeParam=list(mean=mu[1, ], var=Sigma[1, , ]),
                        EdgeParam=list(mean=mu[2, ], var=Sigma[2, , ]))

  res <- list(blockProp=blockProp, connectParam=connectParam, emissionParam=emissionParam)
  return(res)
}

#-------------------------------------------------------------------------------
# E step of the VEM algorithm
#-------------------------------------------------------------------------------
veStep <- function(scoreMat, blockProp, connectParam, emissionParam, tauOld, directed, tauTol=1e-4, etaTol=tauTol){

  # scoreMat <- mat_S; directed <- FALSE; tauTol <- etaTol <- 1e-4; tauOld <- initQuanti$tau[[blockNb]]
  # blockProp <- thetaHat$blockProp; connectParam <- thetaHat$connectParam; emissionParam <- thetaHat$emissionParam

  # Dimensions
  N <- nrow(scoreMat); blockNb <- length(blockProp); n <- nbPairs2n(N, symmetric=!directed)

  # eta
  logPhi <- matrix(0, N, 2)
  logPhi[, 1] <- dmvnorm(scoreMat,
                         mean=emissionParam$noEdgeParam$mean, sigma=emissionParam$noEdgeParam$var)
  logPhi[, 2] <- dmvnorm(scoreMat,
                         mean=emissionParam$EdgeParam$mean, sigma=emissionParam$EdgeParam$var)
  eta <- array(dim=c(N, blockNb, blockNb))
  sapply(1:nBlocks, function(k){sapply(1:nBlocks, function(l){ # k <- 1; l <- 2
    etaTmp <- logPhi * (rep(1, N)%o%c(log(1-connectParam[k, l]), log(connectParam[k, l])))
    etaTmp <- etaTmp - apply(etaTmp, 1, max)
    etaTmp <- exp(etaTmp); etaTmp <- etaTmp / rowSums(etaTmp)
    etaTmp <- etaTmp + etaTol; etaTmp <- etaTmp / rowSums(etaTmp)
    eta[, k, l] <<- etaTmp[, 2]
  })})

  # tau
  logA <- array(dim=c(N, blockNb, blockNb))
  sapply(1:nBlocks, function(k){sapply(1:nBlocks, function(l){ # k <- 1; l <- 2
    logA[, k, l] <<- (1-eta[, k, l])*(log(1-connectParam[k, l]) + logPhi[, 1]) +
      eta[, k, l]*(log(connectParam[k, l]) + logPhi[, 2])
    })})
  indexList <- indices(n, symmetric=!directed)
  tau <- t(sapply(1:n, function(i){ # i <- 3
    indexListIFirst <- which(indexList[, 1]==i)
    indexListISecond <- which(indexList[, 2]==i)
    sapply(1:blockNb, function(k){ # k <- 1
      log(blockProp[k]) +
        sum(logA[indexListIFirst, k, ]*tauOld[indexList[indexListIFirst, 2], ]) +
        sum(logA[indexListISecond, , k]*tauOld[indexList[indexListISecond, 1], ])
      })
    }))
  tau <- tau - apply(tau, 1, max)
  tau <- exp(tau); tau <- tau / rowSums(tau)
  tau <- tau + tauTol; tau <- tau / rowSums(tau)

  # psi
  psi <- rep(0, N)
  sapply(1:blockNb, function(k){sapply(1:blockNb, function(l){ # k <- 1; l <- 2
    psi <<- psi + eta[, k, l] * tau[indexList[, 1], k] *  tau[indexList[, 2], l]
  })})

  res <- list(eta=eta, tau=tau, psi=psi)
  return(res)
}


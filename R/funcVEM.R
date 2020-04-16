###############################################################################
# M step of the VEM algorithm
###############################################################################
mStepNoisySBM <- function(scoreMat, qDist, directed){

  # scoreMat <- scoreMat; qDist <- qDist; directed <- FALSE

  # Dimensions
  d <- ncol(scoreMat); nbBlocks <- ncol(qDist$tau) #N <- nrow(scoreMat);

  # Proportions
  blockProp <- colMeans(qDist$tau)

  # Connection probabilities
  connectParam <- matrix(0, nbBlocks, nbBlocks)
  sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){
    tauVec <- mat2Vect(qDist$tau[, k] %o% qDist$tau[, l], symmetric = !directed, diag = FALSE)
    connectParam[k, l] <<- tauVec %*% qDist$eta[, k, l] / sum(tauVec)
  })})

  # Emission distributions: mu and Sigma
  mu <- matrix(0, 2, d);
  Sigma <- array(dim = c(2, d, d))
  sapply(1:2, function(g){
    mu[g, ] <<- t(qDist$psi[, g]) %*% scoreMat / sum(qDist$psi[, g])
    Sigma[g, , ] <<- t(scoreMat) %*% diag(qDist$psi[, g]) %*% scoreMat / sum(qDist$psi[, g])
    Sigma[g, , ] <<- Sigma[g, , ] - mu[g, ] %o% mu[g, ]
    Sigma[g, , ] <<- .5 * (Sigma[g, , ] + t(Sigma[g, , ]))
  })
  emissionParam <- list(noEdgeParam = list(mean = mu[1, ], var = Sigma[1, , ]),
                        EdgeParam = list(mean = mu[2, ], var = Sigma[2, , ]))

  res <- list(blockProp = blockProp, connectParam = connectParam, emissionParam = emissionParam)
  return(res)
}

###############################################################################
# VE step of the VEM algorithm
###############################################################################
veStepNoisySBM <- function(scoreMat, theta, tauOld, directed, epsilon_tau=1e-4, epsilon_eta=epsilon_tau){

  # scoreMat <- scoreMat; theta <- thetaHat; directed <- FALSE
  # epsilon_tau <- epsilon_eta <- 1e-4; tauOld <- qDist$tau

  # Dimensions
  nbBlocks <- length(theta$blockProp);
  N <- nrow(scoreMat); n <- nbPairs2n(N, symmetric = !directed)
  indexList <- indices(n, symmetric = !directed)

  # log(phi)
  logPhi <- matrix(0, N, 2)
  logPhi[, 1] <- mvtnorm::dmvnorm(scoreMat,
                         mean = theta$emissionParam$noEdgeParam$mean,
                         sigma = theta$emissionParam$noEdgeParam$var, log = TRUE)
  logPhi[, 2] <- mvtnorm::dmvnorm(scoreMat,
                         mean = theta$emissionParam$EdgeParam$mean,
                         sigma = theta$emissionParam$EdgeParam$var, log = TRUE)

  # eta
  eta <- array(dim = c(N, nbBlocks, nbBlocks))
  invisible(sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){ # k <- 1; l <- 2
    etaTmp <- logPhi + (rep(1, N) %o% c(log(1 - theta$connectParam[k, l]), log(theta$connectParam[k, l])))
    etaTmp <- etaTmp - apply(etaTmp, 1, max)
    etaTmp <- exp(etaTmp); etaTmp <- etaTmp / rowSums(etaTmp)
    # etaTmp <- etaTmp + epsilon_eta; etaTmp <- etaTmp / rowSums(etaTmp)
    eta[, k, l] <<- etaTmp[, 2]
  })}))

  # log(A)
  logA <- array(dim = c(N, nbBlocks, nbBlocks))
  sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){ # k <- 1; l <- 2
    logA[, k, l] <<- (1 - eta[, k, l])*(log(1 - theta$connectParam[k, l]) + logPhi[, 1]) +
      eta[, k, l]*(log(theta$connectParam[k, l]) + logPhi[, 2])
    })})

  # tau
  tau <- t(sapply(1:n, function(i){ # i <- 3
    indexListIFirst <- which(indexList[, 1] == i)
    indexListISecond <- which(indexList[, 2] == i)
    sapply(1:nbBlocks, function(k){ # k <- 1
      log(blockProp[k]) +
        sum(logA[indexListIFirst, k, ] * tauOld[indexList[indexListIFirst, 2], ]) +
        sum(logA[indexListISecond, , k] * tauOld[indexList[indexListISecond, 1], ])
      })
    }))
  tau <- tau - apply(tau, 1, max)
  tau <- exp(tau); tau <- tau / rowSums(tau)
  # tau <- tau + epsilon_tau; tau <- tau / rowSums(tau)

  # psi
  psi <- matrix(0, N, 2)
  sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){ # k <- 1; l <- 2
    psi[, 1] <<- psi[, 1] + eta[, k, l] * tau[indexList[, 1], k] *  tau[indexList[, 2], l]
  })})
  psi[,2]  = 1-psi[,1]

  res <- list(eta = eta, tau = tau, psi = psi, logPhi = logPhi)
  return(res)
}

###############################################################################
# Computation of the lowerbound
###############################################################################
lowerBoundNoisySBM <- function(scoreMat,theta,qDist,directed){

  # scoreMat <- scoreMat; theta <- thetaHat; qDist <- qDist

  # Dimensions
  nbBlocks <- length(theta$blockProp);
  N <- nrow(qDist$eta); n <- nbPairs2n(N, symmetric=!directed)
  indexList <- indices(n, symmetric=!directed)

  # Blocks
  espLogpZ <- sum(qDist$tau%*%log(theta$blockProp))
  HqZ <- -sum(qDist$tau*log(qDist$tau + (qDist$tau==0)))

  # Network
  tauArray <- array(dim=c(N, nbBlocks, nbBlocks))
  sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){
    tauArray[, k, l] <<- qDist$tau[indexList[, 1], k] * qDist$tau[indexList[, 2], l]
  })})
  logConnectParam <- log(theta$connectParam); log1_ConnectParam <- log(1-theta$connectParam)
  espLogpG <- sum(sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){
    sum(tauArray[, k, l] * (
      qDist$eta[, k, l] * logConnectParam[k, l] + (1 - qDist$eta[, k, l]) * log1_ConnectParam[k, l]
      ))
  })}))
  HqG <- sum(sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){
    sum(tauArray[, k, l] * (
      qDist$eta[, k, l] * log(qDist$eta[, k, l] + (qDist$eta[, k, l]==0)) +
        (1 - qDist$eta[, k, l]) * log(1 - qDist$eta[, k, l] + (qDist$eta[, k, l]==1))
      ))
  })}))

  # Scores
  espLogpS <- sum(sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){
    sum(tauArray[, k, l] * (
      qDist$eta[, k, l] * theta$logPhi[, 1] + (1 - qDist$eta[, k, l]) * theta$logPhi[, 2]
      ))
  })}))

  # Entropy and lower bound
  entropy <- HqZ + HqG
  lowerBound <- espLogpZ + espLogpG + espLogpS + entropy
  res <- list(espLogpZ=espLogpZ, HqZ=HqZ, espLogpG=espLogpG, HqG=HqG, espLogpS=espLogpZ,
              entropy=entropy, lowerBound=lowerBound)
  return(res)

}

###############################################################################
#---------- Distance on tau--------------------------------------
###############################################################################

distTau  <- function(tau,tauOld)
{
  return(sqrt(sum(as.vector(tau - tauOld)^2)))
}

# distListTheta = function(theta,thetaOld){
#
#   M <- length(theta)
#   D <- sum(vapply(1:M,function(m){sum((theta[[m]] - thetaOld[[m]])^2)},1))
#   return(D)
#
# }


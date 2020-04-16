###############################################################################
#-------------------------------------------------------------------------------
# M step of the VEM algorithm
#-------------------------------------------------------------------------------
###############################################################################
mStepNoisySBM <- function(scoreMat, psi, tau, eta, directed){

  # scoreMat <- mat_S; psi <- initQuanti$Psi; directed <- FALSE
  # nBlocks <- 3; tau <- initQuanti$tau[[nBlocks]]; eta <- initQuanti$eta0[[nBlocks]];

  # Dimensions
  d <- ncol(scoreMat); nbBlocks <- ncol(tau) #N <- nrow(scoreMat);

  # Proportions
  blockProp <- colMeans(tau)

  # Connection probabilities
  connectParam <- matrix(0, nbBlocks, nbBlocks)
  sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){
    tauVec <- mat2Vect(tau[, k] %o% tau[, l], symmetric = !directed, diag = FALSE)
    connectParam[k, l] <<- tauVec %*% eta[, k, l] / sum(tauVec)
  })})

  # Emission distributions: mu and Sigma
  mu <- matrix(0, 2, d);
  Sigma <- array(dim = c(2, d, d))
  sapply(1:2, function(g){
    mu[g, ] <<- t(psi[, g]) %*% scoreMat / sum(psi[, g])
    Sigma[g, , ] <<- t(scoreMat) %*% diag(psi[, g]) %*% scoreMat / sum(psi[, g])
    Sigma[g, , ] <<- Sigma[g, , ] - mu[g, ] %o% mu[g, ]
    Sigma[g, , ] <<- .5 * (Sigma[g, , ] + t(Sigma[g, , ]))
  })
  emissionParam <- list(noEdgeParam = list(mean = mu[1, ], var = Sigma[1, , ]),
                        EdgeParam = list(mean = mu[2, ], var = Sigma[2, , ]))

  res <- list(blockProp = blockProp, connectParam = connectParam, emissionParam = emissionParam)
  return(res)
}
###############################################################################
#-------------------------------------------------------------------------------
# VE step of the VEM algorithm
#-------------------------------------------------------------------------------
###############################################################################
veStepNoisySBM <- function(scoreMat, blockProp, connectParam, emissionParam, tauOld, directed, tauTol=1e-4, etaTol=tauTol){

  # scoreMat <- mat_S; directed <- FALSE; tauTol <- etaTol <- 1e-4; tauOld <- initQuanti$tau[[blockNb]]
  # blockProp <- thetaHat$blockProp; connectParam <- thetaHat$connectParam; emissionParam <- thetaHat$emissionParam

  # Dimensions
  N <- nrow(scoreMat); nbBlocks <- length(blockProp); n <- nbPairs2n(N, symmetric = !directed)
  indexList <- indices(n, symmetric = !directed)

  # log(phi)
  logPhi <- matrix(0, N, 2)
  logPhi[, 1] <- mvtnorm::dmvnorm(scoreMat,
                         mean = emissionParam$noEdgeParam$mean, sigma = emissionParam$noEdgeParam$var, log=TRUE)
  logPhi[, 2] <- mvtnorm::dmvnorm(scoreMat,
                         mean = emissionParam$EdgeParam$mean, sigma = emissionParam$EdgeParam$var, log = TRUE)

  # eta
  eta <- array(dim = c(N, nbBlocks, nbBlocks))
  sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){ # k <- 1; l <- 2
    etaTmp <- logPhi * (rep(1, N) %o% c(log(1 - connectParam[k, l]), log(connectParam[k, l])))
    etaTmp <- etaTmp - apply(etaTmp, 1, max)
    etaTmp <- exp(etaTmp); etaTmp <- etaTmp / rowSums(etaTmp)
    etaTmp <- etaTmp + etaTol; etaTmp <- etaTmp / rowSums(etaTmp)
    eta[, k, l] <<- etaTmp[, 2]
  })})

  # tau
  logA <- array(dim = c(N, nbBlocks, nbBlocks))
  sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){ # k <- 1; l <- 2
    logA[, k, l] <<- (1 - eta[, k, l])*(log(1 - connectParam[k, l]) + logPhi[, 1]) +
      eta[, k, l]*(log(connectParam[k, l]) + logPhi[, 2])
    })})
  tau <- t(sapply(1:n, function(i){ # i <- 3
    indexListIFirst <- which(indexList[, 1] == i)
    indexListISecond <- which(indexList[, 2] == i)
    sapply(1:nbBlocks, function(k){ # k <- 1
      log(blockProp[k]) +
        sum(logA[indexListIFirst, k, ]*tauOld[indexList[indexListIFirst, 2], ]) +
        sum(logA[indexListISecond, , k]*tauOld[indexList[indexListISecond, 1], ])
      })
    }))
  tau <- tau - apply(tau, 1, max)
  tau <- exp(tau); tau <- tau / rowSums(tau)
  tau <- tau + tauTol; tau <- tau / rowSums(tau)

  # psi
  psi <- matrix(0, N,2)
  sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){ # k <- 1; l <- 2
    psi[,1] <<- psi[,1] + eta[, k, l] * tau[indexList[, 1], k] *  tau[indexList[, 2], l]
  })})
  psi[,2]  = 1-psi[,1]

  res <- list(eta = eta, tau = tau, psi = psi, logPhi = logPhi)
  return(res)
}

###############################################################################
#-------------------------------------------------------------------------------
# Computation of the lowerbound
#-------------------------------------------------------------------------------
###############################################################################
borneInfNoisySBM <- function(scoreMat,theta,tau,etaTol,directed){

  # scoreMat <- mat_S
  # blockProp <- thetaHat$blockProp; connectParam <- thetaHat$connectParam; emissionParam <- thetaHat$emissionParam
  # eta <- qDist$eta; tau <- qDist$tau; logPhi <- qDist$logPhi

  nbBlocks <- length(theta$blockProp);   N <- nrow(eta); n <- nbPairs2n(N, symmetric=!directed)
  indexList <- indices(n, symmetric=!directed)

  emissionParam <- theta$emissionParam
  connectParam  <- theta$connectParam

  #---- mise à jour log Phi
  logPhi <- matrix(0, N, 2)
  logPhi[, 1] <- mvtnorm::dmvnorm(scoreMat,
                                  mean = emissionParam$noEdgeParam$mean, sigma = emissionParam$noEdgeParam$var,log = TRUE)
  logPhi[, 2] <- mvtnorm::dmvnorm(scoreMat,
                                  mean = emissionParam$EdgeParam$mean, sigma = emissionParam$EdgeParam$var,log = TRUE)

  #---- mise à jour eta
  eta <- array(dim = c(N, nbBlocks, nbBlocks))
  sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){ # k <- 1; l <- 2
    etaTmp <- logPhi * (rep(1, N) %o% c(log(1 - connectParam[k, l]), log(connectParam[k, l])))
    etaTmp <- etaTmp - apply(etaTmp, 1, max)
    etaTmp <- exp(etaTmp); etaTmp <- etaTmp / rowSums(etaTmp)
    etaTmp <- etaTmp + etaTol; etaTmp <- etaTmp / rowSums(etaTmp)
    eta[, k, l] <<- etaTmp[, 2]
  })})

  #----- mise à jour A

 logA <- array(dim = c(N, nbBlocks, nbBlocks))
  sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){ # k <- 1; l <- 2
    logA[, k, l] <<- (1 - eta[, k, l])*(log(1 - connectParam[k, l]) + logPhi[, 1]) +
      eta[, k, l]*(log(connectParam[k, l]) + logPhi[, 2])
  })})

  #---- calcul Borne Inf
  coeff_sym <- 0.5*(1 - directed) + 1 * (directed)

  D1 <-  sum(tau %*% log(theta$blockProp))  - sum(tau * log(tau))
  D2 <- coeff_sym*(sum(sapply(1:nbBlocks, function(k){ sapply(1:nbBlocks, function(l){
    t(tau[,k]) %*% as.matrix(vect2Mat(logA[,k,l],!directed)) %*% tau[,l]})
  }))) # emission
  D3 <- coeff_sym*(sum(sapply(1:nbBlocks, function(k){ sapply(1:nbBlocks, function(l){
    t(tau[,k]) %*% as.matrix(vect2Mat(eta[, k, l]*log(eta[, k, l]) +  (1 - eta[, k, l])*log(1 - eta[, k, l]),!directed)) %*%
      tau[,l]})})))

  borne_inf <- D1 +  D2 - D3

  return(borne_inf)
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


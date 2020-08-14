veStepScoreSBM <- function(scoreMat, theta, tauOld, directed, epsilon_tau=1e-4, epsilon_eta=epsilon_tau){

  # theta <- thetaHat; epsilon_tau=1e-6; epsilon_eta=epsilon_tau; tauOld <- initDist$tau[[K]]

  # Dimensions
  N <- nrow(scoreMat); nbBlocks <- length(theta$blockProp); n <- nbPairs2n(N, symmetric = !directed)
  indexList <- indices(n = n, symmetric = !directed)

  # log(phi)
  logPhi <- cbind(dmvnorm(scoreMat, mean = theta$emissionParam$noEdgeParam$mean,
                        sigma = theta$emissionParam$noEdgeParam$var, log = TRUE),
                  dmvnorm(scoreMat, mean = theta$emissionParam$EdgeParam$mean,
                        sigma = theta$emissionParam$EdgeParam$var, log = TRUE))

  # eta & log(A)
  logA <- eta <- array(dim=c(N, nbBlocks, nbBlocks))
  sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){
    logAtmp <- cbind(logPhi[, 1] + log(1 - theta$connectParam[k, l]),
                     logPhi[, 2] + log(theta$connectParam[k, l]))
    etaTmp <- logAtmp - apply(logAtmp, 1, max)
    etaTmp <- exp(etaTmp); etaTmp <- etaTmp / rowSums(etaTmp)
    # etaTmp <- etaTmp + epsilon_eta; etaTmp <- etaTmp / rowSums(etaTmp)
    logA[, k, l] <<- etaTmp[, 1] * logAtmp[, 1] + etaTmp[, 2] * logAtmp[, 2]
    eta[, k, l] <<- etaTmp[, 2]
  })})

  # tau
  tau <- t(sapply(1:n, function(i){
    indexIfirst <<- which(indexList[, 1]==i);
    indexIsecond <<- which(indexList[, 2]==i)
    sapply(1:nbBlocks, function(k){
      log(theta$blockProp[k]) +
        sum(logA[indexIfirst, k, ] * tauOld[indexList[indexIfirst, 2], ]) +
        sum(logA[indexIsecond, k, ] * tauOld[indexList[indexIsecond, 1], ])
      })
    }))
  tau <- tau - apply(tau, 1, max)
  tau <- exp(tau); tau <- tau / rowSums(tau)
  # tau <- tau + epsilon_tau; tau <- tau / rowSums(tau)

  # psi
  psi <- rep(0, N)
  sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){
    psi <<- psi + tau[indexList[, 1], k] * tau[indexList[, 2], l] * eta[, k, l]
    })})
  psi <- cbind(psi, 1-psi)

  res <- list(eta = eta, tau = tau, psi = psi, logPhi = logPhi, logA=logA)
  return(res)
}

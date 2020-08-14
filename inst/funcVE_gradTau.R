###############################################################################
# VE step of the VEM algorithm
# TO KEEP: gradient formulas
###############################################################################
veStepScoreSBM <- function(scoreMat, theta, qDistOld, directed, estimOptions = list()){

  # theta <- thetaHat; directed <- FALSE
  # epsilon_tau <- epsilon_eta <- 1e-4; tauOld <- qDist$tau
  currentOptions <- list(
    maxIterVE = 100 ,
    valStopCrit = 1e-6,
    tauTol = 2 * .Machine$double.eps,
    etaTol = 2 * .Machine$double.eps
  )
  currentOptions[names(estimOptions)] <- estimOptions
  qDist <- qDistOld
  lb <- lowerBoundScoreSBM(scoreMat,theta,qDist,directed)
  print(paste('   veOld', lb$lowerBound))
  # print(unlist(lb))

  noConvergence = 0
  # Dimensions
  nbBlocks <- length(theta$blockProp);
  N <- nrow(scoreMat); n <- nbPairs2n(N, symmetric = !directed)
  indexList <- indices(n, symmetric = !directed)

  #-----------------------------------------------------------------------------
  # log(phi)
  logPhi <- matrix(0, N, 2)
  logPhi[, 1] <- mvtnorm::dmvnorm(scoreMat,
                         mean = theta$emissionParam$noEdgeParam$mean,
                         sigma = theta$emissionParam$noEdgeParam$var, log = TRUE)
  logPhi[, 2] <- mvtnorm::dmvnorm(scoreMat,
                         mean = theta$emissionParam$EdgeParam$mean,
                         sigma = theta$emissionParam$EdgeParam$var, log = TRUE)

  #-----------------------------------------------------------------------------
  # eta
  eta <- array(dim = c(N, nbBlocks, nbBlocks))
  invisible(sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){ # k <- 1; l <- 2
    etaTmp <- logPhi + (rep(1, N) %o% c(log(1 - theta$connectParam[k, l] + (theta$connectParam[k, l] == 1)),
                                        log(theta$connectParam[k, l] + (theta$connectParam[k,l] == 0))))
    etaTmp <- etaTmp - apply(etaTmp, 1, max)
    etaTmp <- exp(etaTmp); etaTmp <- etaTmp / rowSums(etaTmp)
    etaTmp <- etaTmp + currentOptions$etaTol; etaTmp <- etaTmp / rowSums(etaTmp)
    eta[, k, l] <<- etaTmp[, 2]
  })}))
  lbOld <- lb
  qDist$eta <- eta; lb <- lowerBoundScoreSBM(scoreMat,theta,qDist,directed)
  print(paste('   eta  ', lb$lowerBound))
  # print(unlist(lb))
  if(lb$lowerBound < lbOld$lowerBound){browser()}

  #-----------------------------------------------------------------------------
  # log(A)
  logA <- array(dim = c(N, nbBlocks, nbBlocks))
  sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){ # k <- 1; l <- 2
    logA[, k, l] <<- (1 - eta[, k, l])*(log(1 - theta$connectParam[k, l] + (theta$connectParam[k,l] == 1)) + logPhi[, 1] - log(1-eta[,k,l])) +
      eta[, k, l]*(log(theta$connectParam[k, l] + (theta$connectParam[k,l] == 0)) + logPhi[, 2] -  log(eta[,k,l]))
    })})

  #-----------------------------------------------------------------------------
  # tau (fixed point)
  lowerBoundTau <- function(tau){
    lowerBoundTau <- sum(tau%*%log(theta$blockProp))
    for(ij in 1:dim(logA)[1]){
      lowerBoundTau <- lowerBoundTau + sum(logA[ij, , ] * (tau[indexList[ij, 1], ]%o%tau[indexList[ij, 2], ]))
      }
    lowerBoundTau <- lowerBoundTau - sum(tau * log(tau + (tau==0)))
    return(lowerBoundTau)
    }
  tau <- qDistOld$tau
  if (nbBlocks > 1){
    iterVE <- 0;  stopVE <- 0
    while ((iterVE < currentOptions$maxIterVE) & (stopVE == 0)) {

      iterVE <- iterVE + 1
      tauOld <- tau;
      tau <- t(sapply(1:n, function(i){ # i <- 3
        indexListIFirst <- which(indexList[, 1] == i)
        indexListISecond <- which(indexList[, 2] == i)
        sapply(1:nbBlocks, function(k){ # k <- 1
          log(theta$blockProp[k]) +
            sum(logA[indexListIFirst, k, ] * tauOld[indexList[indexListIFirst, 2], ]) +
            sum(logA[indexListISecond, , k] * tauOld[indexList[indexListISecond, 1], ])
        })
      }))
      tau <- tau - apply(tau, 1, max)
      tau <- exp(tau); tau <- tau / rowSums(tau)
      tau <- tau + currentOptions$tauTol; tau <- tau / rowSums(tau)

      dTau <- distTau(tau,tauOld)
      if (dTau < currentOptions$valStopCrit)   {stopVE <- 1}
      #print(c(iterVE,dTau))
      if (iterVE == currentOptions$maxIterVE) {noConvergence = noConvergence + 1}
    }
  }
  if(lowerBoundTau(tau) < lowerBoundTau(tauOld)){tau <- tauOld}

  # #-----------------------------------------------------------------------------
  # # tau (gradient descent)
  # negLowerBoundTau <- function(tau1vec){
  #   tau1 <- matrix(tau1vec, n, (nbBlocks-1))
  #   tau <- cbind(tau1, 1-rowSums(tau1))
  #   lowerBoundTau <- sum(tau%*%log(theta$blockProp))
  #   for(ij in 1:dim(logA)[1]){
  #     lowerBoundTau <- lowerBoundTau + sum(logA[ij, , ] * (tau[indexList[ij, 1], ]%o%tau[indexList[ij, 2], ]))
  #     }
  #   lowerBoundTau <- lowerBoundTau - sum(tau * log(tau + (tau==0)))
  #   return(-lowerBoundTau)
  #   }
  # negLowerBoundTauGrad <- function(tau1vec){
  #   tau1 <- matrix(tau1vec, n, (nbBlocks-1))
  #   tau <- cbind(tau1, 1-rowSums(tau1))
  #   grad <- t(sapply(1:n, function(i){ # i <- 3
  #     indexListIFirst <- which(indexList[, 1] == i)
  #     indexListISecond <- which(indexList[, 2] == i)
  #     sapply(1:nbBlocks, function(k){ # k <- 1
  #       log(theta$blockProp[k]) +
  #         sum(logA[indexListIFirst, k, ] * tau[indexList[indexListIFirst, 2], ]) +
  #         sum(logA[indexListISecond, , k] * tau[indexList[indexListISecond, 1], ])
  #       })
  #     }))
  #   grad <- grad - 1 - log(tau)
  #   grad <- grad[, 1:(nbBlocks-1)] - grad[, nbBlocks]%o%rep(1, (nbBlocks-1))
  #   return(-as.vector(grad))
  # }
  # tau1vec <- as.vector(tauOld[, -nbBlocks])
  # tau1vec <- optim(tau1vec, fn=negLowerBoundTau, gr=negLowerBoundTauGrad, method='L-BFGS-B',
  #       lower=rep(currentOptions$tauTol, length(tau1vec)), upper=rep(1-currentOptions$tauTol, length(tau1vec)))$par
  # tau1 <- matrix(tau1vec, n, (nbBlocks-1))

  # #-----------------------------------------------------------------------------
  # # tau (gradient descent : log(tau))
  # negLowerBoundLogTau <- function(logTau1vec){
  #   tau1 <- matrix(exp(logTau1vec), n, (nbBlocks-1))
  #   tau <- cbind(tau1, 1-rowSums(tau1))
  #   lowerBoundTau <- sum(tau%*%log(theta$blockProp))
  #   for(ij in 1:dim(logA)[1]){
  #     lowerBoundTau <- lowerBoundTau + sum(logA[ij, , ] * (tau[indexList[ij, 1], ]%o%tau[indexList[ij, 2], ]))
  #   }
  #   lowerBoundTau <- lowerBoundTau - sum(tau * log(tau + (tau==0)))
  #   return(-lowerBoundTau)
  # }
  # negLowerBoundLogTauGrad <- function(logTau1vec){
  #   tau1 <- matrix(exp(logTau1vec), n, (nbBlocks-1))
  #   tau <- cbind(tau1, 1-rowSums(tau1))
  #   grad <- t(sapply(1:n, function(i){ # i <- 3
  #     indexListIFirst <- which(indexList[, 1] == i)
  #     indexListISecond <- which(indexList[, 2] == i)
  #     sapply(1:nbBlocks, function(k){ # k <- 1
  #       log(theta$blockProp[k]) +
  #         sum(logA[indexListIFirst, k, ] * tau[indexList[indexListIFirst, 2], ]) +
  #         sum(logA[indexListISecond, , k] * tau[indexList[indexListISecond, 1], ])
  #     })
  #   }))
  #   grad <- grad - 1 - log(tau)
  #   grad <- grad[, 1:(nbBlocks-1)] - grad[, nbBlocks]%o%rep(1, (nbBlocks-1))
  #   grad <- grad / tau1
  #   return(-as.vector(grad))
  # }
  # logTau1vec <- as.vector(log(tauOld[, -nbBlocks]))
  # logTau1vec <- optim(logTau1vec, fn=negLowerBoundLogTau, gr=negLowerBoundLogTauGrad, method='L-BFGS-B',
  #                     upper=rep(0, length(logTau1vec)), control=list(maxit=1))$par
  # tau1 <- matrix(exp(logTau1vec), n, (nbBlocks-1))
  #
  # # Completing and smoothing tau's
  # tau <- cbind(tau1, 1-rowSums(tau1))
  # tau <- tau + currentOptions$tauTol
  # tau <- tau/rowSums(tau)

  qDist$tau <- tau
  print(paste('   tau  ', lowerBoundScoreSBM(scoreMat,theta,qDist,directed)$lowerBound))

  #-----------------------------------------------------------------------------
  # psi
  psi <- matrix(0, N, 2)
  sapply(1:nbBlocks, function(k){sapply(1:nbBlocks, function(l){ # k <- 1; l <- 2
    psi[, 2] <<- psi[, 2] + eta[, k, l] * tau[indexList[, 1], k] *  tau[indexList[, 2], l]
  })})
  psi[, 1]  = 1 - psi[, 2]

  qDist$psi <- psi
  print(paste('   psi  ', lowerBoundScoreSBM(scoreMat,theta,qDist,directed)$lowerBound))

  qDist <- list(eta = eta, tau = tau, psi = psi)
  return(qDist)
}


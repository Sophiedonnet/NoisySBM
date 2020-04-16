rm(list=ls())

library(NoisySBM)
library(tensor)
# source('D:/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Stephane_Robin/NOISY_Papier_Package/Package/NoisySBM/R/funcVEM.R', encoding = 'UTF-8')
# source('D:/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Stephane_Robin/NOISY_Papier_Package/Package/NoisySBM/R/tools.R', encoding = 'UTF-8')
# source('D:/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Stephane_Robin/NOISY_Papier_Package/Package/NoisySBM/R/tools.R', encoding = 'UTF-8')
source('./R/funcVEM.R', encoding = 'UTF-8')
source('./R/tools.R', encoding = 'UTF-8')
library(mclust)

# Parms
nNodes  <- 20 #60
blockProp <- c(1/3,1/2,1/6)
nBlocks   <- length(blockProp) # SR: 'mixtureParam' -> 'blockProp'
connectParam <- matrix(rbeta(nBlocks^2,1.5,1.5 ),nBlocks,nBlocks)
connectParam <- 0.5*(connectParam + t(connectParam))
emissionParam <- list()
d <- 4
emissionParam$noEdgeParam <- list(mean = rep(0,d),var = diag(0.1,nrow = d,ncol = d))
emissionParam$EdgeParam <- list( mean = 1:d,var =  diag(0.1,nrow = d,ncol = d))

# Data
directed <- FALSE
dataSim <- rNoisySBM(nNodes, directed = TRUE, blockProp, connectParam, emissionParam, seed = NULL)

# Init
K <- 3
scoreList <- dataSim$noisyNetworks
initDist <- initInferenceNoisySBM(scoreList, directed)
scoreMat <- scoreList2scoreMat(scoreList, symmetric=!directed)
thetaInit <- mStepNoisySBM(scoreMat=scoreMat,
                          qDist=list(psi=initDist$psi, tau=initDist$tau[[K]], eta=initDist$eta[[K]]),
                          directed=directed)

# VEM
thetaHat <- thetaInit
# maxIterVE <- NULL; epsilon_tau <- 1e-4; epsilon_eta <- 2 * .Machine$double.eps
maxIterVEM <- 100; J <- rep(0, 2*maxIterVEM)
for(iter in 1:maxIterVEM){
  qDist <- veStepNoisySBM(scoreMat=scoreMat, theta=thetaHat, tauOld=initDist$tau[[K]], directed=directed)
  crit <- lowerBoundNoisySBM(scoreMat=scoreMat,theta=thetaHat,qDist=qDist,directed)
  J[2*iter-1] <- crit$lowerBound
  thetaHat <- mStepNoisySBM(scoreMat=scoreMat, qDist=qDist, directed=directed)
  crit <- lowerBoundNoisySBM(scoreMat=scoreMat,theta=thetaHat,qDist=qDist,directed)
  J[2*iter] <- crit$lowerBound
  plot(J[1:iter], col=rep(1:2, iter), pch=20, type='b')
}
print(unlist(crit))

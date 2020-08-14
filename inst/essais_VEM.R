library(ScoreSBM)
nbNodes  = 60
directed = TRUE
blockProp = c(1/3,1/2,1/6)
nbBlocks   = length(blockProp)
connectParam <- matrix(rbeta(nbBlocks^2,1.5,1.5 ),nbBlocks,nbBlocks)
connectParam <- 0.5*(connectParam + t(connectParam))
emissionParam <- list()
nbScores <- 4
emissionParam$noEdgeParam = list(mean = rep(0,nbScores),var = diag(0.1,nrow = nbScores,ncol = nbScores))
emissionParam$edgeParam = list( mean = 1:nbScores,var =  diag(0.1,nrow = nbScores,ncol = nbScores))
data1 <- rScoreSBM(nbNodes,directed, blockProp,connectParam,emissionParam,seed = NULL)


N <- nbNodes*(nbNodes - 1)*(0.5 + 0.5*directed)
scoreList <- data1$scoreNetworks
source('R/tools.R')
scoreMat <- sapply(1:nbScores , function(q) {mat2Vect(scoreList[[q]], symmetric = !directed, diag = F)})

library(mclust)
initAll <- initInferenceScoreSBM(scoreList, directed)
nbModels <- length(initAll$tau)
nbBlocksAll = vapply(1:nbModels,function(m){ncol(initAll$tau[[m]])},1)

rangeK <- length(initAll$tau)
bestK <- which.max(initAll$ICL)
bestK <- 3
psi <- initAll$psi
tauBestK <- initAll$tau[[bestK]]
etaBestK <- initAll$eta[[bestK]]

qDist = list(tau = tauBestK,eta = etaBestK,psi = initAll$psi)
# essai etape M


source('R/funcVEM.R')
theta <- mStepScoreSBM(scoreMat, qDist , directed)
J1 <- lowerBoundScoreSBM(scoreMat,theta,qDist,directed)
J1
L1 <- J1$lowerBound
# essai etap VE (1 iter)

qDist <- veStepScoreSBM(scoreMat, theta,qDist$tau, directed)
L2 <- lowerBoundScoreSBM(scoreMat,theta,qDist,directed)$lowerBound

print(c(L1,L2))
##################################################################
#--------------------- VE -----------------------------
##################################################################
initBestK <- list(psi = initAll$psi)
initBestK$tau = initAll$tau[[bestK]]
initBestK$eta = initAll$eta[[bestK]]
init = initBestK
resVEM <- VEMScoreSBM(scoreMat, directed, init,monitoring = list(lowerBound = TRUE),estimOptions = list(verbosity = 0,maxIterVE = 100,maxIterVEM = 100))

L <- lowerBoundScoreSBM(scoreMat,resVEM$theta,resVEM$qDist,directed)

plot(resVEM$lowerBound,type = 'l')

##################################################################
#--------------------- ESTIMATE ALL  -----------------------------
##################################################################
res_ALL <- estimateScoreSBM(scoreList,directed)

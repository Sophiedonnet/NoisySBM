
nbNodes  = 100
directed = TRUE
blockProp = c(1/3,1/2,1/6)
nbBlocks   = length(blockProp)
connectParam <- matrix(rbeta(nbBlocks^2,1.5,1.5 ),nbBlocks,nbBlocks)
connectParam <- 0.5*(connectParam + t(connectParam))
emissionParam <- list()
nbScores <- 2
emissionParam$noEdgeParam = list(mean = rep(0,nbScores),var = diag(0.1,nrow = nbScores,ncol = nbScores))
emissionParam$EdgeParam = list( mean = 1:nbScores,var =  diag(0.1,nrow = nbScores,ncol = nbScores))
data1 <- rNoisySBM(nbNodes,directed, blockProp,connectParam,emissionParam,seed = NULL)


N <- nbNodes*(nbNodes - 1)*(0.5 + 0.5*directed)
scoreList <- data1$noisyNetworks
scoreMat <- sapply(1:nbScores , function(q) {mat2Vect(scoreList[[q]], symmetric = !directed, diag = F)})


initAll <- initInferenceNoisySBM(scoreList, directed)
rangeK <- length(initAll$tau)
bestK <- which.max(initAll$ICL)
psi <- initAll$psi
tauBestK <- initAll$tau[[bestK]]
etaBestK <- initAll$eta[[bestK]]

qDist = list(tau=tauBestK,eta=etaBestK,psi = initAll$psi)
# essai etape M
theta <- mStepNoisySBM(scoreMat, qDist , directed)
J1 <- lowerBoundNoisySBM(scoreMat,theta,qDist,directed)
# essai etap VE (1 iter)
tauTol <- 0.001
etaTol <- 0.001
maxIterVE <- 100
qDist <- veStepNoisySBM(scoreMat, theta,qDist$tau, directed, tauTol, etaTol,maxIterVE,valStopCrit = 0.00001)
J2 <- lowerBoundNoisySBM(scoreMat,theta,qDist,directed)


##################################################################
#--------------------- VE -----------------------------
##################################################################
initBestK <- list(psi = initAll$psi)
initBestK$tau = initAll$tau[[bestK]]
initBestK$eta = initAll$eta[[bestK]]
init = initBestK
resVEM <- VEMNoisySBM(scoreMat, directed, init,monitoring = list(lowerBound = TRUE),maxIterVE = 100 ,
                      maxIterVEM = 10)


plot(resVEM$lowerBound,type = 'l')




##################################################################
#--------------------- BorneInf -----------------------------
##################################################################
initBestK <- list(psi = initAll$psi)
initBestK$tau = initAll$tau[[bestK]]
initBestK$eta = initAll$eta[[bestK]]
init = initBestK
resVEM <- VEMNoisySBM(scoreMat, directed, init,monitoring = list(lowerBound = TRUE),maxIterVE = 100 ,
                      maxIterVEM = 1000)


plot(resVEM$lowerBound,type = 'l')

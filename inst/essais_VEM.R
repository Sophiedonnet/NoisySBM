
nbNodes  = 60
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

# essai etape M
resM <- mStepNoisySBM(scoreMat, initAll$psi, tauBestK , etaBestK , directed)

# essai etap VE (1 iter)
tauTol <- 0.001
etaTol <- 0.001
tau0  <- tauBestK
resVE <- veStepNoisySBM(scoreMat, resM$blockProp, resM$connectParam, resM$emissionParam, tau0, directed, tauTol = tauTol, etaTol = etaTol)
J <- borneInfNoisySBM(scoreMat,resM,resVE$tau,etaTol,directed)



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

#



##################################################################
#--------------------- BorneInf -----------------------------
##################################################################
initBestK <- list(psi = initAll$psi)
initBestK$tau = initAll$tau[[bestK]]
initBestK$eta = initAll$eta[[bestK]]
init = initBestK
resVEM <- VEMNoisySBM(scoreMat, directed, init,monitoring = list(lowerBound = TRUE),maxIterVE = 100 ,
                      maxIterVEM = 10)


plot(resVEM$lowerBound,type = 'l')

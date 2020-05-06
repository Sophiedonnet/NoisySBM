context("testing VEM steps")

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




##################################################################
#---------------------INITIALISATION -----------------------------
##################################################################
initAll <- initInferenceNoisySBM(scoreList, directed)
rangeK <- length(initAll$tau)
bestK <- which.max(initAll$ICL)
psi <- initAll$psi
tauBestK <- initAll$tau[[bestK]]
etaBestK <- initAll$eta[[bestK]]

test_that("The initialisation is giving  elements of right size", {

  ########### several models :  must be a list
  expect_type(initAll,"list")

  ########### Dimension and type of the objects
  expect_true(is.matrix(initAll$psi))
  expect_equal(nrow(initAll$psi),N)
  expect_equal(ncol(initAll$psi),2)

  expect_type(initAll$tau,"list")
  expect_type(initAll$eta,"list")
  expect_equal(length(initAll$eta),rangeK)

  ######### Check the size of each matrix


  expect_true(is.array(etaBestK))
  expect_equal(dim(etaBestK),c(N,bestK,bestK))

  expect_true(is.matrix(tauBestK))
  expect_equal(nrow(tauBestK),nbNodes)
  expect_equal(ncol(tauBestK),bestK)
}
)


##################################################################
#---------------------stepM -----------------------------
##################################################################
qDist <- list(tau = tauBestK,eta = etaBestK,psi = initAll$psi)
resM <- mStepNoisySBM(scoreMat, qDist, directed)

test_that("M step is working as espected", {


  expect_true(is.vector(resM$blockProp))
  expect_equal(c(length(resM$blockProp),sum(resM$blockProp)),c(bestK,1))


  expect_true(is.matrix(resM$connectParam))
  expect_equal(dim(resM$connectParam),c(bestK,bestK))

  expect_type(resM$emissionParam,"list")
  expect_equal(length(resM$emissionParam),2)

  expect_equal(length(resM$emissionParam$EdgeParam),2)
  expect_equal(length(resM$emissionParam$noEdgeParam),2)
  expect_equal(dim(resM$emissionParam$noEdgeParam$var),c(nbScores,nbScores))
  expect_equal(dim(resM$emissionParam$EdgeParam$var),c(nbScores,nbScores))

  expect_equal(sum(diag(resM$emissionParam$EdgeParam$var) <= 0),0)
  expect_equal(sum(diag(resM$emissionParam$noEdgeParam$var) <= 0),0)
  }
)
##################################################################
#---------------------stepVE -----------------------------
##################################################################

resVE <- veStepNoisySBM(scoreMat, theta = resM,qDist$tau, directed)
test_that("VE step is working as espected", {

  expect_type(resVE,"list")

  expect_true(is.matrix(resVE$psi))
  expect_equal(dim(resVE$psi),c(N,2))

  expect_true(is.array(resVE$eta))
  expect_equal(dim(resVE$eta),c(N,bestK,bestK))

  expect_true(is.matrix(resVE$tau))
  expect_equal(dim(resVE$tau),c(nbNodes,bestK))
}
)

##################################################################
#--------------------- BorneInf -----------------------------
##################################################################
J <- lowerBoundNoisySBM(scoreMat,theta = resM,qDist,directed)$lowerBound
test_that("Lower Bound is working well", {

   expect_true(is.numeric(J))
 }
)

##################################################################
#--------------------- VEM-----------------------------
##################################################################
resVEM <-  VEMNoisySBM(scoreMat, directed, qDist,monitoring = list(lowerBound = TRUE),estimOptions = list(verbosity = 0,maxIterVE = 10,maxIterVEM = 10))


test_that("Lower Bound is working well", {

  expect_true(is.numeric(J))
}
)

context("testing sampling rNoisySBM")

nbNodes  = 100
directed = TRUE
blockProp = c(1/3,1/2,1/6)
nbBlocks   = length(blockProp)
connectParam <- matrix(rbeta(nbBlocks^2,1.5,1.5 ),nbBlocks,nbBlocks)
connectParam <- 0.5*(connectParam + t(connectParam))
emissionParam <- list()
nbScores <- 4
emissionParam$noEdgeParam = list(mean = rep(0,nbScores),var = diag(0.1,nrow = nbScores,ncol = nbScores))
emissionParam$edgeParam = list( mean = 1:nbScores,var =  diag(0.1,nrow = nbScores,ncol = nbScores))
data1 <- rNoisySBM(nbNodes,directed = TRUE, blockProp,connectParam,emissionParam,seed = NULL)


S <- data1$noisyNetworks
G <- data1$trueNetwork
Z <- data1$memberships

test_that("The simulated object of correct type and  dimensions", {

  ########### dataSim must be a list
  expect_type(data1,"list")

  ########### S must be a list
  expect_type(S,"list")
  expect_equal(length(S),nbScores)

  ######### Check the size of each matrix
  test_dim <- 1;
  for (net in 1:nbScores) {
    dim_net <- dim(S[[net]])
    dim_theo <- c(nbNodes,nbNodes)
    test_dim <- test_dim * prod(dim_net == dim_theo)
  }
  expect_equal(test_dim,1)

  ######### Check the symmetry
  test_symmetry = 1

  if (directed) {
    test_symmetry <- test_symmetry *  isSymmetric(G) ;
    for (net in 1:nbScores) {
      test_symmetry <- test_symmetry * isSymmetric(S[[net]])
    }
  }
  expect_equal(test_symmetry,1)

  ######### Check the type of the content of each matrix
  test_G <- as.numeric(all( G %in% c(0,1)))
  expect_equal(test_G,1)
}
)







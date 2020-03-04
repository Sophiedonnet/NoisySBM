context("testing sampling rNoisySBM")

nNodes  = 100
symmetric = TRUE
mixtureParam = c(1/3,1/2,1/6)
nBlocks   = length(mixtureParam)
connectParam <- matrix(rbeta(nBlocks^2,1.5,1.5 ),nBlocks,nBlocks)
connectParam <- 0.5*(connectParam + t(connectParam))
emissionParam <- list()
d <- 4
emissionParam$noEdgeParam = list(mean=rep(0,d),var = diag(0.1,nrow = d,ncol = d))
emissionParam$EdgeParam = list( mean= 1:d,var =  diag(0.1,nrow = d,ncol = d))
data1 <- rNoisySBM(nNodes, symmetric = TRUE, mixtureParam,connectParam,emissionParam,seed = NULL)

S <- data1$noisyNetworks
G <- data1$trueNetwork
Z <- data1$memberships

test_that("The simulated object of correct type and  dimensions", {

  ########### dataSim must be a list
  expect_type(data1,"list")

  ########### S must be a list
  expect_type(S,"list")
  expect_equal(length(S),d)

  ######### Check the size of each matrix
  test_dim <- 1;
  for (net in 1:d) {
    dim_net <- dim(S[[net]])
    dim_theo <- c(nNodes,nNodes)
    test_dim <- test_dim * prod(dim_net == dim_theo)
  }
  expect_equal(test_dim,1)

  ######### Check the symmetry
  test_symmetry = 1

  if (symmetric) {
    test_symmetry <- test_symmetry *  isSymmetric(G) ;
    for (net in 1:d) {
      test_symmetry <- test_symmetry * isSymmetric(S[[net]])
    }
  }
  expect_equal(test_symmetry,1)

  ######### Check the type of the content of each matrix
  test_G <- as.numeric(all( G %in% c(0,1)))
  expect_equal(test_G,1)
}
)








#------------------------------------------------------------
#--------------------Intialisation of the inference procedure
#------------------------------------------------------------
initInference <- function(mat_S,symmetric){

  N <- nrow(mat_S)
  #--------  Mclust on the d scores dyad per dyad
  param_gm <- mclust::Mclust(mat_S, G = 2, verbose = FALSE)
  Psi <- param_gm$z
  G <- param_gm$classification - 1
  mu  <- param_gm$parameters$mean
  test_G <- rowMeans(t(mu)) #identify G = 0  and G =1
  if (test_G[1] > test_G[2]) {
    Psi <- Psi[,c(2,1)]
    G = 1 - G
  }


  #------------------ init of SBM parameters
  membership_type <- ifelse(symmetric, "SBM_sym", "SBM")
  param_sbm <- blockmodels::BM_bernoulli(membership_type, adj = vect2Mat(G, symmetric),
                                       plotting = '',
                                       verbosity = 0)
  param_sbm$estimate()
  Kmax <- length(param_sbm$memberships)
  tau_init <-  lapply(1:Kmax, function(K){param_sbm$memberships[[K]]$Z})
  eta0 <-  lapply(1:Kmax,  function(K){array(rep(Psi[, 1], K * K), c(N, K, K))})
  #------------------------------------------------------------

  res <- list(Psi = Psi,tau  = tau_init,eta0  = eta0, ICL  = param_sbm$ICL, G=G)
  return(res)

}



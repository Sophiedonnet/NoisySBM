VEM_NoisySBM <- function(data, symmetric, K,
    maxIterVE = NULL ,
    maxIterVEM = NULL,
    explorFact = 1.5,
    epsilon_tau = 1e-4,
    epsilon_eta = 2 * .Machine$double.eps) {


    #------------------------------------------------------------
    #eps <- 2 * .Machine$double.eps
    valStopCrit <- 1e-6
    #------------------------- Data under various format
    S <- data  # a list at that point.
    d <- length(S);
    n <- nrow(S[[1]])

    array_S <- array(dim = c(d, n, n))
    invisible(sapply(1:d, function(h) {array_S[h, ,] <<- S[[h]]}))
    mat_S <- sapply(1:d, function(q) {mat2Vect(S[[q]], symmetric = symmetric, diag = F)})

    N <- nrow(mat_S)
    transfo_indices <- indices(n, symmetric)
    A = array(0, dim = c(N, K, K))


    #--------------------Intialisation of parameters eta and G (G a deux classes arete ou pas)

    param_gm <- mclust::Mclust(mat_S, G = 2, verbose = FALSE)
    Psi <- param_gm$z
    G <- param_gm$classification - 1
    mu  <- param_gm$parameters$mean
    #identify G = 0  and G =1
    test_G <- rowMeans(t(mu))
    if (test_G[1] > test_G[2]) {
      Psi <- Psi[,c(2,1)]
      G = 1 - G
    }
    eta0 <- array(rep(Psi[, 1], K * K), c(N, K, K)) #ok
    eta1 <- 1 - eta0
    #------------------ init of SBM parameters
    membership_type = ifelse(symmetric, "SBM_sym", "SBM")
    param_sbm <- blockmodels::BM_bernoulli(membership_type, adj = vect2Mat(G, symmetric),
        plotting = '',
        verbosity = 0,
        ncores = 1,
        exploration_factor = explorFact
      )
    param_sbm$estimate()
    tau_init <- tau_hat <- param_sbm$memberships[[K]]$Z
    theta <- list()
    theta$gamma <- param_sbm$model_parameters[[K]]$pi
    theta$pi <- colMeans(tau_hat)

    #------------------ init of emission parameters
    theta$mu0 <- mu[,1]
    theta$mu1 <- mu[,2]
    theta$var0 <- cov(mat_S[G == 0,])
    theta$var1 <- cov(mat_S[G == 1,])



    #------------- nb of max  iteration
    if (is.null(maxIterVE)) { maxIterVE = 1000}
    if (is.null(maxIterVEM)){maxIterVEM  =  1000}
    stopCrit <- 0
    iterVEM <- 0

    J  <- numeric(3 * maxIterVEM)


    thetaOld <- theta

    ######################
    # Algo begins
    #####################
    noConvergence = 0
    while (iterVEM < maxIterVEM & stopCrit == 0)
    {
      iterVEM <- iterVEM + 1




      #--------------------------------  VE step
      #---------------------- COMPUTING the Eta and eta bar = psi
      EA <- computeEta0A(mat_S,theta,epsilon_eta)
      eta0 <- EA$eta0
      A <- EA$A
      eta1 <- 1 - eta0;

      #---- calcul des tau
      iterVE <- 0;  stopVE <- 0
      tau_old <-  tau_hat
      tau_new <- matrix(0, n, K)


      while ((iterVE < maxIterVE) & (stopVE == 0)) {

        #----------------- for each i = 1 ... n
        invisible(sapply(1:n, function(i) {
          # m1.i = ensemble des paires dont i est le premier element
          m1.i = which(transfo_indices[, 1] == i)
          if (length(m1.i) > 0) { B1.i = array(A[m1.i, ,], dim = c(length(m1.i), K, K))}
          # m2.i = ensemble des paires dont i est le second element
          m2.i = which(transfo_indices[, 2] == i)
          # Transpoition de A quand i est le second element de la paire m
          if (length(m2.i) > 0) {
            B2.i = array(A[m2.i, ,], dim = c(length(m2.i), K, K))
            sapply(1:length(m2.i), function(m) {
              B2.i[m, ,] <<- t(B2.i[m, ,])
            })
          }
          B.i = array(dim = c(length(m1.i) + length(m2.i), K, K))
          if (length(m1.i) > 0) { B.i[1:length(m1.i), ,] = B1.i  }
          if (length(m2.i) > 0) { B.i[(length(m1.i) + 1):(length(m1.i) + length(m2.i)), ,] = B2.i}
          #  Calcul de tau
          tau.j = tau_old[-i,]
          tau.i = rep(0, K)
          invisible(sapply(1:K, function(k) {tau.i[k] <<- log(theta$pi[k]) + sum(tau.j * B.i[, k,])}))
          tau.i = tau.i - max(tau.i)
          tau.i[which(tau.i < -100)] = -100
          tau.i = exp(tau.i)
          tau.i = tau.i / sum(tau.i)
          tau.i = tau.i + 1e-4
          tau.i = tau.i / sum(tau.i)
          tau_new[i,] <<- tau.i
        })) #% ------------ end for sapply
      iterVE <- iterVE + 1
      if (distTau(tau_new,tau_old) < valStopCrit)   {stopVE <- 1}
      if (iterVE == maxIterVE) {noConvergence = noConvergence + 1}
      tau_old = tau_new
      }
      tau_hat = tau_new
      J[(3 * iterVEM) - 2] =  borneInfV2(mat_S,theta,tau_hat,epsilon_eta,symmetric)


      #--------------------------------   M step

      #  Calcul gamma:
      Gamma_hat_num <- matrix(nrow = K, ncol = K)
      Gamma_hat_num <-
        sapply(1:K, function(k) {
          sapply(1:K, function(l) {
            t(tau_hat[, k]) %*% vect2Mat(eta1[, k, l], symmetric) %*% tau_hat[, l]
          })
        })
      Gamma_hat_denum <- (t(tau_hat) %*% (matrix(1, n, n) - diag(1, n)) %*% tau_hat)
      theta$gamma <- Gamma_hat_num / Gamma_hat_denum
      #- Calcul pi
      theta$pi <- colMeans(tau_hat)

      J[(3 * iterVEM) - 1] =  borneInfV2(mat_S,theta,tau_hat,epsilon_eta,symmetric)

      #  Calcul de phi les parametres des densités
      Psi0_hat =  rep(0, N)
      invisible(sapply(1:(n - 1), function(j) {
        sapply((j + 1):n, function(i) {
          ij = which((transfo_indices[, 1] == i) & (transfo_indices[, 2] == j))
          # cat(i, j, ij, transfo_indices[ij, ], '\n')
          Psi0_hat[ij] <<-
            (t(tau_hat[i, ]) %*% eta0[ij, ,] %*% tau_hat[j, ])[1, 1]
          #Psi1_hat[ij] <<-  (t(tau_hat[i,]) %*% eta1[ij, , ] %*% tau_hat[j,])[1, 1]
        })
      }))
      #calcul de phi les parametres des densités
      Psi1_hat = 1 - Psi0_hat

      #   Calcul de  paramEmission
      theta$mu1 <-  as.vector(Psi1_hat %*% mat_S / sum(Psi1_hat))
      theta$mu0 <-  as.vector(Psi0_hat %*% mat_S / sum(Psi0_hat))
      theta$var1 <- t(mat_S) %*% diag(Psi1_hat) %*% mat_S / sum(Psi1_hat) - theta$mu1 %o% theta$mu1
      theta$var0 <- t(mat_S) %*% diag(Psi0_hat) %*% mat_S / sum(Psi0_hat) - theta$mu0 %o% theta$mu0


      J[(3 * iterVEM) ] =  borneInfV2(mat_S,theta,tau_hat,epsilon_eta,symmetric)

      if (distListTheta(theta, thetaOld) < valStopCrit){ stopCrit <- 1}
      thetaOld <- theta

plot(J[J != 0],type = 'l')
    }
    ############" reorder
    ord <- order(diag(theta$gamma), decreasing  = TRUE)
    theta$pi <- theta$pi[ord]
    theta$gamma <- theta$gamma[ord,ord]

    J  = J[J != 0]
    plot(J,type='l')
    output <- list(tau_init  = tau_init[,ord], tau  = tau_hat[,ord],borne_inf = J)
    output$theta <- theta
    output$Psi0 <- Psi0_hat
    output$Psi1 <- Psi1_hat
    return(output)




  }

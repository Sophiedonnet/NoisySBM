############################################################
############ Fonctions de bases à utiliser dans R ####################
############ Ex : "Main_function_last" ####################
############ Manipulation de vecteurs, etc ####################
############################################################

#############################################################
#Construit une matrice symétrique à partir d'un vecteur
#############################################################
vect_to_mat <- function(V, symmetric,diag=F)
{

  N <- length(V)

  if (symmetric == TRUE){
    n <- ifelse(diag, (-1 + sqrt(1 + 8*N))/2 , (1 + sqrt(1 + 8 * N))/2)
    M <- matrix(0,n,n)
    M[lower.tri(M, diag = diag)] <- V
    B <- t(M)
    diag(B) <- 0
    M <- B + M
  }else{
    n <- ifelse(diag, sqrt(N) , (1 + sqrt(1 + 4*N))/2)
    M <- matrix(0,n,n)
    where <- (lower.tri(M, diag = diag) +  upper.tri(M, diag = diag)) > 0
    M[where] <- V
  }
  return(M)
}

###########################################################################
#Construit un vecteur à partir d'une matrice triangulaire inférieure
###########################################################################
mat_to_vect <- function(M, symmetric, diag=F)
{
  if (symmetric) {
    where <- lower.tri(M, diag = diag)
  } else {
    where <- (lower.tri(M, diag=diag) +  upper.tri(M, diag = diag)) > 0
  }
  V <- M[where]
  return(V)
}

############################################################################
#Passage d'une matrice binaire à un vecteur des indices et inversement
############################################################################
mat_bin_to_vect_ind <- function(Z)
{
  K <- ncol(Z);
  Z_ind <- Z %*% (1:K)
  return(Z_ind)
}

vect_ind_to_mat_bin <- function(Zind)
{
  K <- length(unique(Zind))
  n <- length(Zind)
  Zind = sample(1:4,n,replace=T)
  M <-vapply(1:K,function(k){Mk <- rep(0,n);  Mk[which(Zind==k)] = 1; return(Mk)},rep(0,n))
  return(M)
}

###########################################################################
#recuperer les indices de la matrice triangulaire inferieure (dans le cas symmetric) , ou hors diag si besoin
###########################################################################

indices <- function(n, symmetric,diag=F)
{
  if (symmetric){
    N  <- (n * (n - 1)/2)*(diag==F) + (n * (n + 1)/2)*(diag==T)
  }else{
    N <- (n * n) * (diag == T) + (n * n - n) * (diag == F)
  }
  S <- vect_to_mat(c(1:N), symmetric,diag = diag)

  if (symmetric) {S[upper.tri(S)] = 0}
  res <- which(S!=0,arr.ind=T)
  return(res)
}



#############################################################################
#############################################################################
# Indice_mat <- function(n,V)
# {
#   i <- (V-1)%/%n +1
#   j <- V%%n + n * ((V%%n)==0)
#   return(cbind(matrix(c(i,j),ncol=2)))
# }


#############################################################################
# #############################################################################
# fun_test<-function(a_i_b_i){
#   return(exp(seq(log(a_i_b_i["a"]),log(a_i_b_i["b"]),length.out = a_i_b_i["k"])))
# }
#############################################################################
#############################################################################
##########################################################################
# Calcul de la borne inf
#### non vérifées
borne_inf <- function(tau_hat,Pi_hat,A,eta0,eta1,symmetric){
  coeff_sym <- 0.5*(symmetric) + 1 * (1-symmetric)
  borne_inf <- sum(tau_hat %*% log(Pi_hat))  - sum(tau_hat * log(tau_hat)) +
    coeff_sym*(sum(sapply(1:K, function(k){ sapply(1:K, function(l){
      t(tau_hat[,k]) %*% as.matrix(vect_to_mat(A[,k,l],symmetric)) %*% tau_hat[,l]})
    }))) -
    coeff_sym*(sum(sapply(1:K, function(k){ sapply(1:K, function(l){
      t(tau_hat[,k]) %*%
        as.matrix(vect_to_mat(eta0[, k, l]*log(eta0[, k, l]) +
            eta1[, k, l]*log(eta1[, k, l]),symmetric)) %*%
        tau_hat[,l]})
    })))
  return(borne_inf)
}


#### non vérifées
borneInfV2 <- function(mat_S,theta,tau_hat,epsilon_eta,symmetric){

  EA <- computeEta0A(mat_S,theta,epsilon_eta)
  eta0 <- EA$eta0
  A <- EA$A
  eta1 <- 1 - eta0;


  coeff_sym <- 0.5*(symmetric) + 1 * (1-symmetric)
  D1 <-  sum(tau_hat %*% log(theta$pi))  - sum(tau_hat * log(tau_hat))
  D2 <- coeff_sym*(sum(sapply(1:K, function(k){ sapply(1:K, function(l){
    t(tau_hat[,k]) %*% as.matrix(vect_to_mat(A[,k,l],symmetric)) %*% tau_hat[,l]})
  })))
  D3 <- coeff_sym*(sum(sapply(1:K, function(k){ sapply(1:K, function(l){
    t(tau_hat[,k]) %*% as.matrix(vect_to_mat(eta0[, k, l]*log(eta0[, k, l]) +  eta1[, k, l]*log(eta1[, k, l]),symmetric)) %*%
      tau_hat[,l]})})))

  borne_inf <- D1 +  D2 - D3

  return(borne_inf)
}

#### non vérifées
computeEta0A <- function(mat_S,theta,epsilon_eta){
  lfu <- cbind(
    mvtnorm::dmvnorm(mat_S, theta$mu0, theta$var0, log = TRUE),
    mvtnorm::dmvnorm(mat_S, theta$mu1, theta$var1, log = TRUE)
  )
  dlfu = lfu[, 1] - lfu[, 2]
  dlfu[which(abs(dlfu) > 100)] = sign(dlfu[which(abs(dlfu) > 100)]) * 100 # Troncature pour eviter les pb numeriques
  eta1 = 1 / (1 + exp(dlfu) %o% ((1 - theta$gamma) / theta$gamma))
  eta0 = 1 - eta1
  # Lissage des proba a conditionnelles
  eta0 = eta0 + epsilon_eta ; eta1 = eta1 + epsilon_eta ; eta = eta0 + eta1;
  eta0 = eta0 / eta;
  A <- eta0 * (rep(1, N) %o% log(1 - theta$gamma)) + eta0 * lfu[, 1] + eta1 * (rep(1, N) %o% log(theta$gamma)) + eta1 * lfu[, 2]

  return(list(eta0 = eta0,A=A))

}

#--------------------------
#----------------------  distance on tau

distTau  <- function(tau,tauOld)
{
  Q <- length(tau)
  vdis <- sapply(1:Q,function(q){
    return(sqrt(sum(as.vector(tau[[q]] - tauOld[[q]])^2)))
  })
  return(sum(vdis))
}

distListTheta = function(theta,thetaOld){

  M <- length(theta)
  D <- sum(vapply(1:M,function(m){sum((theta[[m]] - thetaOld[[m]])^2)},1))
  return(D)

}


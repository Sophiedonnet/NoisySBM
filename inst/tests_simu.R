nNodes  = 60
directed = TRUE
mixtureParam = c(1/3,1/2,1/6)
nBlocks   = length(mixtureParam)
connectParam <- matrix(rbeta(nBlocks^2,1.5,1.5 ),nBlocks,nBlocks)
connectParam <- 0.5*(connectParam + t(connectParam))
emissionParam <- list()
d <- 4
emissionParam$noEdgeParam = list(mean=rep(0,d),var = diag(0.1,nrow = d,ncol = d))
emissionParam$EdgeParam = list( mean= 1:d,var =  diag(0.1,nrow = d,ncol = d))
dataSim <- rNoisySBM(nNodes,symmetric = TRUE, mixtureParam,connectParam,emissionParam,seed = NULL)


symmetric = TRUE
maxIterVE = NULL ;
maxIterVEM = NULL;
explorFact = 1.5;
epsilon_tau = 1e-4;
epsilon_eta = 2 * .Machine$double.eps
K  = 3
data <- dataSim$noisyNetworks



source('D:/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Stephane_Robin/NOISY_Papier_Package/Package/NoisySBM/R/VEM_NoisySBM.R', encoding = 'UTF-8')
source('D:/WORK_ALL/RECHERCHE/TRAVAUX_RECHERCHE/Stephane_Robin/NOISY_Papier_Package/Package/NoisySBM/R/tools.R', encoding = 'UTF-8')
library(mclust)
res_VEM <- VEM_NoisySBM(data, symmetric, K,
                                    maxIterVE = NULL ,
                                    maxIterVEM = NULL,
                                    explorFact = 1.5,
                                    epsilon_tau = 1e-4,
                                    epsilon_eta = 2 * .Machine$double.eps)

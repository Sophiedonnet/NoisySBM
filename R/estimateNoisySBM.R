#' Inference of the Noisy SBM  Model
#'
#' \code{estimateNoisySBM} performs the estimation and model selection for the NoisySBM model.
#' @param scoreList   :  list of Scores for each dyad  of an underlying network
#' @param directed    :  if true the inference network is directed. Default value = FALSE.
#' @param estimOptions : tunes the optimization process (see details below)
#' @param monitoring : specifies if the lowerBound along the VEM iterations is saved (monitoring = list(lowerBound = TRUE))
#' @details  The list of parameters \code{estimOptions} essentially tunes the optimization process and the variational EM algorithm, with the following parameters
#'  \itemize{
#'  \item{"verbosity"}{ controls the verbosity of the procedure (0 or 1). Default is 1.}
#'  \item{"exploreFactor"}{ controls the exploration of the number of groups in the initialization step, Default is 1.5}
#'  \item{"nbBlocksRange"}{ minimal and maximal number or blocks explored. Default is c(1,Inf)}
#'  \item{"nbCores"}{ integer for number of cores used. Default is 1. }
#'  \item{"maxIterVE"}{ Maximum number of iterations in the VE step. Default is 100.}
#'  \item{"maxIterVE"}{ Maximum number of iterations in the VEM algorithm. Default is 1000.}
#'  \item{"tauTol"}{ Tolerance in the VEM algorithm. Default value \code{tauTol = 2 * .Machine$double.eps}}
#'  \item{"etaTol"}{ Tolerance in the VEM algorithm. Default value \code{etaTol = 2 * .Machine$double.eps}}
#'  \item{"valCritStop"}{ Algorithm stops when the difference on successive tau is below valCritStop. Default value is  1e-6}
#'  }

#' @return The output is a list of estimated models, each one corresponding to a number of blocks.
#' @details Each element of the output list contains the following quantites:
#'  \itemize{
#'  \item{"theta"}{ estimated parameters}
#'  \item{"nbBlocks"}{ number of blocks in the underlying network}
#'  \item{"qDist"}{ variational approximation of the the conditional distribution q(G,Z | data)}
#'  \item{"ICL": }{Integrated Likelihood Criterion}
#'  \item{"pen": }{Penalty term for model selection (included in the ICL)}
#'  \item{"lowerBound": }{Lowerbound along the VEM iterations (provided if \code{monitoring$lowerBound = TRUE})}
#'  }
#' @examples
#' nbNodes  <- 60
#' directed <- FALSE
#' blockProp <- c(1/3,1/2,1/6)
#' nbBlocks   <- length(blockProp)
#' connectParam <- matrix(rbeta(nbBlocks^2,1.5,1.5 ),nbBlocks,nbBlocks)
#' connectParam <- 0.5*(connectParam + t(connectParam))
#' emissionParam <- list()
#' nbScores <- 4
#' emissionParam$noEdgeParam <- list(mean=rep(0,nbScores));
#' emissionParam$noEdgeParam$var <- diag(0.1,nrow = nbScores,ncol = nbScores)
#' emissionParam$EdgeParam <- list( mean= 1:nbScores)
#' emissionParam$EdgeParam$var <-  diag(0.1,nrow = nbScores,ncol = nbScores)
#' dataSim <- rNoisySBM(nbNodes,directed = TRUE, blockProp,connectParam,emissionParam,seed = NULL)
#' scoreList <- dataSim$noisyNetworks
#' resEstim <- estimateNoisySBM(scoreList,directed)
#' @importFrom parallel mclapply
#' @importFrom pbmcapply pbmclapply
#' @export


estimateNoisySBM = function(scoreList,directed = FALSE, estimOptions=list(), monitoring = list()){

  currentOptions <- list(
    verbosity     = 1,
    explorFactor  = 1.5,
    nbBlocksRange = c(1,Inf),
    nbCores       = 1,
    maxIterVE = 100 ,
    maxIterVEM = 1000,
    tauTol = 2 * .Machine$double.eps,
    valCritStop = 1e-6,
    etaTol = 2 * .Machine$double.eps
  )
  currentMonitoring = list(
    lowerBound = TRUE
  )

  ## Current options are default expect for those passed by the user
  currentOptions[names(estimOptions)] <- estimOptions
  currentMonitoring[names(monitoring)] <- monitoring


  #------------------------ transform Data
  nbScores <- length(scoreList)
  scoreMat <- sapply(1:nbScores , function(q) {mat2Vect(scoreList[[q]], symmetric = !directed, diag = F)})

  #--------------------------------- initialisation .
  if (currentOptions$verbosity > 0) { print("-------------- Initialization ----------- ")}
  initAll <- initInferenceNoisySBM(scoreList, directed,currentOptions)

  #-------------------------------------------------------------------------
  #--------------------ESTIMATION ALL MODELS -------------------------------
  #------------------------------------------------------------------------
  nbModels <- length(initAll$tau)
  nbBlocksAll = vapply(1:nbModels,function(m){ncol(initAll$tau[[m]])},1)
  indModels <- (1:nbModels)[nbBlocksAll <= currentOptions$nbBlocksRange[2]]
  indModels <- (1:nbModels)[nbBlocksAll >= currentOptions$nbBlocksRange[1]]
  nbBlocksAll <- nbBlocksAll[indModels]
  if (currentOptions$verbosity > 0) { print("-------------- Estimation of the models ----------- ")}

  init_m  <- list(psi = initAll$psi)
  Estim_m <- function(m){
    if ((currentOptions$verbosity > 0 ) & (currentOptions$nbCores == 1)) {
      print(paste("-------------- Estimation for",m ,"blocks---- ----- ",sep = ' '))
    }
    init_m$tau = initAll$tau[[m]]
    init_m$eta = initAll$eta[[m]]
    resVEM_m <- VEMNoisySBM(scoreMat, directed, init_m,currentOptions,currentMonitoring)
    if ((currentOptions$verbosity > 0 ) & (currentOptions$nbCores == 1)) {
      print(paste("ICL = ",resVEM_m$ICL,sep = ' '))
    }

    return(resVEM_m)
  }

  listRes <- pbmclapply(indModels,Estim_m,mc.cores = currentOptions$nbCores);

  #-------------------------------------------------------------------------
  #--------------------REORDERING  MODELS -------------------------------
  #------------------------------------------------------------------------
  ordModels  <- order(sapply(indModels,function(m){listRes[[m]]$ICL}),decreasing  = TRUE)
  res <- lapply(ordModels, function(m){listRes[[m]]})
  return(res)

}








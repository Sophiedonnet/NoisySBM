#' Inference of the Noisy SBM  Model
#
#' @param data        :  list of d Noisy versions of a network
#' @param directed    :  if true the inference network is directed. Default value = FALSE.
#'
#' @details The function \code{multipartiteBM} selects the better numbers of blocks in each FG (with a penalized likelihood criterion). The model selection is performed with a forward backward strategy and the likelihood of each model is maximized with a variational EM).
#' @return a list of estimated parameters for the different models ordered by decreasing ICL. If save \code{= FALSE}, the length is of length 1 (only the better model is returned).
#' \describe{




#'   }
#'   }
#'}
#' @examples
#' nNodes <- 100
#' @importFrom parallel mclapply
#' @export


estimateNoisySBM = function(scoreList,directed = FALSE, estimOptions=list(), monitoring = list()){

  currentOptions <- list(
    verbosity     = 0,
    plot          = character(0),
    explorFactor  = 1.5,
    nbBlocksRange = c(1,Inf),
    nbCores       = 1,
    maxIterVE = 100 ,
    maxIterVEM = 1000,
    tauTol = 2 * .Machine$double.eps,
    valCritStop = 1e-6,
    etaTol = 2 * .Machine$double.eps
  )
  currentMonitoring = list(lowerBound = FALSE)

  ## Current options are default expect for those passed by the user
  currentOptions[names(estimOptions)] <- estimOptions
  currentMonitoring[names(monitoring)] <- monitoring


  #------------------------ transform Data
  scoreMat <- sapply(1:nbScores , function(q) {mat2Vect(scoreList[[q]], symmetric = !directed, diag = F)})

  #--------------------------------- initialisation .
  if (currentOptions$verbosity > 0) { print("-------------- Initialization ----------- ")}
  initAll <- initInferenceNoisySBM(scoreList, directed,currentOptions)

  #-------------------------------------------------
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
    return(resVEM_m)
  }

  listRes <- mclapply(indModels,Estim_m,mc.cores = currentOptions$nbCores);

  resModelSelect <- as.data.frame(t(sapply(indModels,function(m){
    u <- c(listRes[[m]]$nbBlocks,listRes[[m]]$ICL, listRes[[m]]$pen)
  })))
  names(resModelSelect) <- c('nbBlocks','ICL','penalty')

}








#' Inference of the Noisy SBM  Model
#
#' @param data        :  list of d Noisy versions of a network
#' @param symmetric   :  if true the inference network is symmetric (non directed),
#' @param nBlocksMin  :  minimum number of Blocks
#' @param nBlocksMax  :  maximum number of Blocks
#' @param nBlocksInit :  initialisation number of Blocks
#' @param save        :  an optional boolean. If TRUE  save the estimated parameters for intermediate visited models. Otherwise, only the better model (in ICL sense) is the ouput. Default value \code{= FALSE}.
#' @param verbose     :  an optional boolean. If  TRUE, display the current step of the search algorithm
#' @param nCores      :  an optional integer specifying the number or cores used for the estimation. Not parallelized on windows. If \code{ncores = NULL}, then half of the cores are used.
#' @param maxIterVE   :  an optional integer  specifying the maximum number of iterations in the VE step of the VEM algorithm. If NULL then default value  \code{= 1000}
#' @param maxIterVEM  :  an optional integer  specifying the maximum number of iterations of the VEM algorithm. If NULL then default value Default value  \code{= 1000}
#' @details The function \code{multipartiteBM} selects the better numbers of blocks in each FG (with a penalized likelihood criterion). The model selection is performed with a forward backward strategy and the likelihood of each model is maximized with a variational EM).
#' @return a list of estimated parameters for the different models ordered by decreasing ICL. If save \code{= FALSE}, the length is of length 1 (only the better model is returned).
#' \describe{
#'   \item{\code{fittedModel}}{contains the results of the inference. \code{res$fittedModel[[1]]}  is a list with fields
#'   \describe{
#'   \item{\code{paramEstim}}{a MBMfit object.}
#'   \item{\code{ICL}}{the penalized likelihood criterion ICL.}
#'   \item{\code{LowerBound}}{the sequence of the varational bound of the likelihood through iterations of the VEM.}
#'   \item{\code{monitoring}}{TRUE if the VEM reached convergence.}
#'   }
#'   }
#'}
#' @examples
#' nNodes <- 100

#' @export


estimateNoisySBM = function(data,symmetric,nBlocksMin = NULL, nBlocksMax = NULL ,nBlocksInit  = NULL ,save=FALSE , verbose = TRUE, nCores = NULL, maxIterVE = NULL , maxIterVEM = NULL) {



  #------------------- Messages  ----------
  if (verbose) {
    nNodes <- nrow(data[[1]]);
    print("------------Nb of Nodes   ------------")
    print(nNodes)
  }
  #------------------------  Bounds of the number of clusters
  if (is.null(nBlocksMin)) {nBlocksMin = 1; print("The minimum number of blocks has been set to 1")}
  if (is.null(nBlocksMax)) {nBlocksMax = 1; print("The maximum number of blocks has been set to 1")}
  #------------------------ nBlocksInit
  if (is.null(nBlocksInit)) { nBlocksInit  = 1 }

  #-------------------------




}

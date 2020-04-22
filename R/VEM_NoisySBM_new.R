VEMNoisySBM <- function(scoreMat, directed, init,
    maxIterVE = NULL ,
    maxIterVEM = NULL,
    explorFact = 1.5,
    tauTol = 1e-4,
    etaTol = 2 * .Machine$double.eps, monitoring = list(lowerBound = FALSE)) {


    #------------------------------------------------------------
    #eps <- 2 * .Machine$double.eps
    if (is.null(maxIterVE))  maxIterVE = 100
    if (is.null(maxIterVEM))  maxIterVEM  =  1000
    valStopCrit <- 1e-6
    stopCrit <- 0
    iterVEM <- 0
    if (monitoring$lowerBound) J  <- numeric(2 * maxIterVEM)

    ############################################
    #------------ initialisation
    ############################################



    ############################################
    #--------------   Algo begins
    ###########################################
    #noConvergence = 0

    qDist = init
    while (iterVEM < maxIterVEM & stopCrit == 0)
    {
      print(iterVEM)
      iterVEM <- iterVEM + 1
      tauCurrent <- qDist$tau
      #---------------------------------------
      #------------  M step ------------------
      #---------------------------------------

      #-------
      theta <- mStepNoisySBM(scoreMat, qDist, directed)
      #------

      if (monitoring$lowerBound) J[(2 * iterVEM) - 1] =   lowerBoundNoisySBM(scoreMat,theta,qDist,directed)$lowerBound



      #-------------------------------------
      #-------------- VE step ----------------
      #---------------------------------------
      qDist <- veStepNoisySBM(scoreMat, theta,tauOld = qDist$tau, directed, tauTol, etaTol,maxIterVE,valStopCrit)
      if (monitoring$lowerBound) J[(2 * iterVEM)] =  lowerBoundNoisySBM(scoreMat,theta,qDist,directed)$lowerBound


      deltaTau <- distTau(tauCurrent, qDist$tau)
      if ( deltaTau < valStopCrit) { stopCrit <- 1}
      print(deltaTau)


    }
    ############################################
    #--------------  End of the algorithm
    ###########################################


    ############################################
    #------- Reorder
    #########################################
    ord <- order(diag(theta$connectParam), decreasing  = TRUE)
    theta$blockProp <- theta$blockProp[ord]
    theta$connectParam <- theta$connectParam[ord,ord]

    output <- list(tau  = tau[,ord], theta = theta)
    if (monitoring$lowerBound) output$lowerBound <- J[1:(2*iterVEM)]

    return(output)




  }

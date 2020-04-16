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

    psi <- init$psi;
    tau <- init$tau;
    eta <- init$eta;

    ############################################
    #--------------   Algo begins
    ###########################################
    noConvergence = 0
    while (iterVEM < maxIterVEM & stopCrit == 0)
    {
      print(iterVEM)
      iterVEM <- iterVEM + 1
      tauCurrent <- tau
      #---------------------------------------
      #------------  M step ------------------
      #---------------------------------------
      resM <- mStepNoisySBM(scoreMat, psi, tau, eta, directed)
      connectParam <- resM$connectParam
      blockProp <- resM$blockProp
      emissionParam <- resM$emissionParam
      theta = list(connectParam = connectParam, blockProp = blockProp, emissionParam = emissionParam)

      if (monitoring$lowerBound) J[(2 * iterVEM) - 1] =  borneInfNoisySBM(scoreMat,theta,tau,etaTol,directed)
      #---------------------------------------
      #-------------- VE step ----------------
      #---------------------------------------
      iterVE <- 0;  stopVE <- 0

      while ((iterVE < maxIterVE) & (stopVE == 0)) {

        iterVE <- iterVE + 1


        tauOld <- tau;
        resVE <- veStepNoisySBM(scoreMat, blockProp, connectParam, emissionParam, tauOld, directed, tauTol = tauTol, etaTol = etaTol)
        tau <- resVE$tau

        dTau <- distTau(tau,tauOld)
        if (dTau < valStopCrit)   {stopVE <- 1}
        print(c(iterVE,dTau))
        if (iterVE == maxIterVE) {noConvergence = noConvergence + 1}
      }

      psi <- resVE$psi
      eta <- resVE$eta
      if (monitoring$lowerBound) J[(2 * iterVEM)] =  borneInfNoisySBM(scoreMat,theta,tau,etaTol,directed)
      #---------------------------------------
      #-------------- VE step ----------------
      #---------------------------------------



      if (distTau(tauCurrent, tau) < valStopCrit) { stopCrit <- 1}



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
    if (monitoring$lowerBound) output$lowerBound <- J

    return(output)




  }

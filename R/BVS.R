
# Internal function that performs one MCMC iteration to sample from the posterior distribution of all unkonws involved in the 
# Bayesian Variable Selection step. For details, see Alexopoulos and Bottolo (2019)


########## BVS ##########

BVS <- function(cors, Y, Bcurr, XB, Xactual, Rcurr, invRcurr, SS, invSS, ZZ, cutoff, 
                response_type, api, bpi, ncat, commoncovs, gammacurr, surecovs, dd, negbinompar, trials)
{
  p_adddel <- 0.9
  p_add <- 0.5
  geom_mean <- 5
  
  ########## comment ##########
  # Hard coding! in the gamma proposal, "p_adddel" is the probability of adding or deleting (0.1 probability of swapping); "p_add" is 
  # the probability of adding (0.5 probability of deleting); "geom_mean" is mean of the geometric distribution mixture component
  
  n <- dim(ZZ)[1]
  m <- dim(ZZ)[2]
  conmu <- convar <- TU <- TL <- matrix(NA, n, m)
  TUprop <- TLprop <- rep(NA, n)
  Ztilde <- ZZ + XB
  tau1 <- 1
  myj <- 0
  myj2 <- 0
  gammaprop <- list()
  for (j in 1 : m)
  {
    indna <- which(!is.na(Y[, j]))
    start <- (j - 1) * (commoncovs + surecovs) + 1
    end <- j * (commoncovs + surecovs)
    mygamma <- which(gammacurr[start : end] > 0)
    
    ########## propose gamma ##########
    
    gammaprop <- Sample_gamma(cors[, j], dim(Xactual)[2] - surecovs, p_adddel, gammacurr[(start + surecovs) : end], p_add, geom_mean)
    loggammapropratio <- gammaprop[[2]]
     
    if (surecovs > 0)
    {
      gammaprop[[1]] <- c(1 : surecovs, gammaprop[[1]] + surecovs)
    }
     
    if (length(gammaprop[[1]] > 0))
    {
      setone <- gammaprop[[1]] + (j - 1) * dim(Xactual)[2]
      setzero <- ((start : end)[-gammaprop[[1]]])
    }else{
      setone <- NULL
      setzero <- (start : end)
    }
    
    if (surecovs > 0)
    {
      gammastar_initial <- gammaprop[[1]][-c(1 : surecovs)]
    }else{
      gammastar_initial <- gammaprop[[1]]
    }
    
    logpriorgammaratio <- lbeta(length(gammastar_initial) + api, (dim(Xactual)[2] - surecovs) - length(gammastar_initial) + bpi) - 
                          lbeta(sum(gammacurr[start : end]) - surecovs + api, dim(Xactual)[2] - surecovs - sum(gammacurr[start : end]) + bpi)
    
    if (length(mygamma) > 0)
    {
      Xgamma <- as.matrix(Xactual[, mygamma])
      Bj <- Bcurr[start : end]
      Bjgamma <- Bj[mygamma]
      XBgamma <- Xgamma %*% Bjgamma
      transXgamma <- t( Xgamma)
    }else{
      Xgamma <- matrix(0, n, 1)
      Bj <- rep(0, commoncovs)
      Bjgamma <-0
      XBgamma <- rep(0, n)
      transXgamma <- t(Xgamma)
    }
    
    if (length(gammaprop[[1]]) > 0)
    {
      if (surecovs > 0)
      {
        gammastar <- c(1 : surecovs, gammastar_initial)
      }else{
        gammastar <- c(gammastar_initial)
      }
      Xgammastar <- as.matrix(Xactual[, gammastar])
      transXgammastar <- t(Xgammastar)
    }else{
      if (surecovs > 0)
      {
        gammastar <- c(1 : surecovs)
        Xgammastar <- as.matrix(Xactual[, gammastar])
        transXgammastar <- t(Xgammastar)
      }else{
        gammastar <- NULL
        Xgammastar <- matrix(0, n, 1)
        transXgammastar <- t(Xgammastar)
      }
    }
    
    ########## propose beta ##########
    
    formu <- Rcurr[j, -j] %*% solve(Rcurr[-j, -j])
    conmu[, j] <- formu %*% t(ZZ[, -j])
    convar[, j] <- Rcurr[j, j] - formu %*% Rcurr[-j, j]
    if (response_type[j] == "Gaussian")
    {
      Ztilde[, j] <- dd[j] * ZZ[, j] + XB[, j]
      const <- invRcurr[j, j] / (dd[j] * dd[j])
      indshelp <- 1 : m
      myTT <- rep(0, n)
      for (kk in indshelp[-j])
      {
        myTT <- myTT + ZZ[, kk] * invRcurr[j, kk] / dd[j]
      }
      if (length(gammastar) > 1) 
      {
        Q <- diag(rep(1 /(tau1), length(gammastar))) + const * transXgammastar %*% Xgammastar
        invQ <- solve(Q)
      }else{
        Q <- as.matrix(1 /(tau1) ) + const * transXgammastar %*% Xgammastar
        invQ <- solve(Q)
      }
      muprop <- invQ %*% transXgammastar %*% (const * Ztilde[, j] + myTT)
      if (length(mygamma) > 1)
      {
        Qold <- diag(rep(1 /(tau1), length(mygamma))) + const * transXgamma %*% Xgamma
        invQold <- solve(Qold)
      }else{
        Qold <- as.matrix(1 /(tau1)) + const * transXgamma %*% Xgamma
        invQold <- solve(as.matrix(1 /(tau1)) + const * transXgamma %*% Xgamma)
      }
      mupropold <- invQold %*% transXgamma %*% (const * Ztilde[, j] + myTT)
      Bprop <- mvrnorm(1, muprop, invQ)
    }else{
      Ztilde[, j] <- ZZ[, j] + XB[, j]
      const <- invRcurr[j, j]
      indshelp <- 1 : m
      myTT <- rep(0, n)
      for (kk in indshelp[-j])
      {
        myTT <- myTT + ZZ[, kk] * invRcurr[j, kk]
      }
      if (length(gammastar) > 1)
      {
        Q <- diag(rep(1 /(tau1), length(gammastar))) + const * transXgammastar %*% Xgammastar
        invQ <- solve(Q)
      }else{
        Q <- as.matrix(1 /(tau1)) + const*transXgammastar%*%Xgammastar
        invQ <- solve(Q)
      }
      muprop <- invQ %*% transXgammastar %*% (const * Ztilde[, j] + myTT)
      if (length(mygamma) > 1)
      {
        Qold <- diag(rep(1 /(tau1), length(mygamma))) + const * transXgamma %*% Xgamma
        invQold <- solve(Qold)
      }else{
        Qold <- as.matrix(1 /(tau1)) + const * transXgamma %*% Xgamma
        invQold <- solve(as.matrix(1 /(tau1))  + const * transXgamma %*% Xgamma)
      }
      mupropold <- invQold %*% transXgamma %*% (const * Ztilde[, j] + myTT)
      Bprop <- mvrnorm(1, muprop, invQ)
    }
    if (length(gammastar) > 0)
    {
      XBprop<-(as.matrix(Xactual[, gammastar]) %*% Bprop)[, 1]
    }else{
      XBprop <- rep(0,n)
    }
        
    ########## update beta and gamma ##########
    
    if (response_type[j] == "Gaussian")
    {
      LikeHolmes <- -0.5 * determinant(Qold, logarithm = TRUE)$modulus + 0.5 * t(mupropold) %*% Qold %*% as.matrix(mupropold) - 0.5 * length(mygamma) * log(tau1)
      LikeHolmesProp <- -0.5 * determinant(Q, logarithm = TRUE)$modulus + 0.5 * t(muprop) %*% Q %*% as.matrix(muprop) - 0.5 * length(gammastar) * log(tau1)
      LogaccRatioCon <- LikeHolmesProp - LikeHolmes + loggammapropratio + logpriorgammaratio
      if (log(runif(1)) < LogaccRatioCon)
      {
        if (length(gammastar) > 0)
        {
          Bcurr[setone] <- Bprop
          gammacurr[setone] <- 1
        }
        gammacurr[setzero] <- 0
        Bcurr[setzero] <- 0
        XBgamma <- XBprop
        XB[, j] <- XBprop
      }
      for (i in 1 : n)
      {
        if (is.na(Y[i, j]))
        {
          eps <- rnorm(1, formu %*% (ZZ[i, -j]), sqrt(Rcurr[j, j] - formu %*% Rcurr[-j, j]))
          ZZ[i,j] <- eps
        }else{
          ZZ[i, j] <- (Y[i, j] - XBgamma[i]) / dd[j]
        }
      }
    }else{
      if (response_type[j] == "binary")
      {
        LikeHolmes <- -0.5 * determinant(Qold, logarithm = TRUE)$modulus + 0.5 * t(mupropold) %*% Qold %*% as.matrix(mupropold) - 0.5*length(mygamma) * log(tau1)
        LikeHolmesProp <- -0.5 * determinant(Q, logarithm = TRUE)$modulus + 0.5 * t(muprop) %*% Q %*% as.matrix(muprop) - 0.5 * length(gammastar) * log(tau1)
        LogaccRatioCon <- LikeHolmesProp - LikeHolmes + loggammapropratio + logpriorgammaratio
        if (log(runif(1)) < LogaccRatioCon)
        {
          if (length(gammastar) > 0)
          {
            Bcurr[setone] <- Bprop
            gammacurr[setone] <- 1
          }
          gammacurr[setzero] <- 0
          Bcurr[setzero] <- 0
          XBgamma <- XBprop
          XB[, j] <- XBprop
        }
        ZZ[, j] <- Ztilde[, j] - XB[, j]
        TU[, j] <- cutoff[[j]][Y[, j] + 2] - XBgamma
        TL[, j] <- cutoff[[j]][Y[, j] + 1] - XBgamma
      }
      
      if (response_type[j] == "ordinal") 
      {
        LikeHolmes <- -0.5 * determinant(Qold, logarithm = TRUE)$modulus + 0.5 * t(mupropold) %*% Qold %*% as.matrix(mupropold) - 0.5 * length(mygamma) * log(tau1)
        LikeHolmesProp <- -0.5 * determinant(Q, logarithm = TRUE)$modulus + 0.5 * t(muprop) %*% Q %*% as.matrix(muprop) - 0.5 * length(gammastar) * log(tau1)
        LogaccRatioCon <- LikeHolmesProp - LikeHolmes + loggammapropratio + logpriorgammaratio
        if (log(runif(1)) < LogaccRatioCon)
        {
          if (length(gammastar) > 0)
          {
            Bcurr[setone] <- Bprop
            gammacurr[setone] <- 1
          }
          gammacurr[setzero] <- 0
          Bcurr[setzero] <- 0
          XBgamma <- XBprop
          XB[, j] <- XBprop
        }
        ZZ[, j] <- Ztilde[, j] - XB[, j]
        TU[, j] <- cutoff[[j]][Y[, j] + 2] - XBgamma
        TL[, j] <- cutoff[[j]][Y[, j] + 1] - XBgamma
      }
      
      if (response_type[j] == "negative binomial")
      {
        myj <- myj+1
        TUprop <- qnorm(pnbinom(Y[, j], negbinompar[myj], 1 - (1 /(1 + exp(-XBprop)))))
        TU[, j] <- qnorm(pnbinom(Y[, j], negbinompar[myj], 1 - (1 /(1 + exp(-XBgamma)))))
        TLprop <- qnorm(pnbinom(Y[, j] - 1, negbinompar[myj], 1 - (1 /(1 + exp(-XBprop)))))
        TL[,j]<- qnorm( pnbinom(Y[, j] - 1, negbinompar[myj], 1 - (1 /(1 + exp(-XBgamma)))))
        yy <- union(which(TU[, j] == TL[, j]), which(TL[, j] == Inf))
        if ((length(yy) > 0))
        {
          TL[, j][yy] <- qnorm(0.9999999999999)
          TU[, j][yy] <- qnorm(0.99999999999999)
          
          ########## comment ##########
          # For some outliers $y_{ik}$, the pdf $F_{k}(y_{ik}|\beta_{k},\theta_{k})$ and $F_{k(y_{ik} - 1|\beta_{k},\theta_{k})$ can
          # become one in R and qnorm(1) = Inf, but TL and TU cannot be equal since they are the truncation points for the Gaussian 
          # distribution from which the latent variables are drawn. Therefore, whenever $F_{k}(y_{ik}|\beta_{k},\theta_{k})$=1 in
          # R, we set $F_{k}(y_{ik}|\beta_{k},\theta_{k})$=0.99999999999999 and $F_{k(y_{ik} - 1|\beta_{k},\theta_{k})$= 0.9999999999999
          
        }
        yy <- union(which(TUprop == TLprop), which(TLprop == Inf))
        if ((length(yy) > 0))
        {
          TLprop[yy] <- qnorm(0.9999999999999)
          TUprop[yy] <- qnorm(0.99999999999999)
          
          ########## comment ##########
          # For some outliers $y_{ik}$, the pdf $F_{k}(y_{ik}|\beta_{k},\theta_{k})$ and $F_{k(y_{ik} - 1|\beta_{k},\theta_{k})$ can
          # become one in R and qnorm(1) = Inf, but TL and TU cannot be equal since they are the truncation points for the Gaussian 
          # distribution from which the latent variables are drawn. Therefore, whenever $F_{k}(y_{ik}|\beta_{k},\theta_{k})$=1 in
          # R, we set $F_{k}(y_{ik}|\beta_{k},\theta_{k})$=0.99999999999999 and $F_{k(y_{ik} - 1|\beta_{k},\theta_{k})$= 0.9999999999999
          
        }
      }
      
      if (response_type[j] == "binomial")
      {
        myj2 <- myj2 + 1
        TUprop <- qnorm(pbinom(Y[, j], trials[myj2], 1 /(1 + exp(-XBprop))))
        TU[, j] <- qnorm(pbinom(Y[, j], trials[myj2], 1 /(1 + exp(-XBgamma))))
        TLprop <- qnorm(pbinom(Y[, j] - 1, trials[myj2], 1 /(1 + exp(-XBprop))))
        TL[, j] <- qnorm(pbinom(Y[, j] - 1, trials[myj2], 1 /(1 + exp(-XBgamma))))
      }
      
      if (response_type[j] == "binomial" || response_type[j] == "negative binomial")
      {
        DiffLikProp <- pnorm((TUprop[indna] - conmu[indna, j]) / sqrt(convar[indna, j])) - pnorm((TLprop[indna] - conmu[indna, j]) / sqrt(convar[indna, j]))
        DiffLikProp[which(DiffLikProp <= 0)] <- pnorm(Inf) - pnorm(8.2)
        
        ########## comment ##########
        # In equation (12) of the main paper (Alexopoulos and Bottolo, 2019), there is the difference of two standard Gaussian cdfs, 
        # but this difference can be come very small (zero in R) if the cdfs are both evaluated at large values. Instead, we compute 
        # pnorm(Inf) - pnorm(8.2) which is the smallest difference that can be computed in R
        
        logLikProp <- sum(log(DiffLikProp))
        DiffLik <- ((pnorm((TU[indna, j] - conmu[indna, j]) / sqrt(convar[indna, j])) - pnorm((TL[indna, j] - conmu[indna, j]) / sqrt(convar[indna, j]))))
        DiffLik[which(DiffLik <= 0)] <- pnorm(Inf) - pnorm(8.2)
        
        ########## comment ##########
        # In equation (12) of the main paper (Alexopoulos and Bottolo, 2019), there is the difference of two standard Gaussian cdfs, 
        # but this difference can be come very small (zero in R) if the cdfs are both evaluated at large values. Instead, we compute 
        # pnorm(Inf) - pnorm(8.2) which is the smallest difference that can be computed in R
        
        logLikOld <- sum(log(DiffLik))
        logPrior <- sum(dnorm(Bprop, 0, sqrt(tau1), log = TRUE)) - sum(dnorm(Bjgamma, 0, sqrt(tau1), log = TRUE))
        logAccRatio <- logPrior + dmvnorm(Bjgamma, mupropold, invQold, log = TRUE) - dmvnorm(Bprop, muprop, invQ, log = TRUE) + logpriorgammaratio
        LogaccRatioDisc <- logAccRatio + logLikProp - logLikOld + loggammapropratio
        if (log(runif(1)) < LogaccRatioDisc)
        {
          if (length(gammastar) > 0)
          {
            Bcurr[setone] <- Bprop
            gammacurr[setone] <- 1
            gammacurr[setzero] <- 0
          }
          Bcurr[setzero] <- 0
          XBgamma <- XBprop
          XB[, j] <- XBprop
          TU[, j] <- TUprop
          TL[, j] <- TLprop
        }
      }
      
      ########## update cut-offs ##########
      
      if (response_type[j] == "ordinal")
      {
        mygamma <- which(gammacurr[start : end] > 0)
        if (length(mygamma) > 0)
        {
          Xgamma <- as.matrix(Xactual[, mygamma])
          Bj <- Bcurr[start : end]
          Bjgamma <- Bj[mygamma]
          XBgamma <- Xgamma %*% Bjgamma
          transXgamma <- t(Xgamma)
        }else{
          Xgamma <- matrix(0, n, 1)
          Bj <- rep(0, commoncovs)
          Bjgamma <- 0
          XBgamma <- rep(0, n)
          transXgamma <- t(Xgamma)
        }
        a <- rep(NA, length(cutoff[[j]]) - 2)
        a[1] <- 0
        a[2 : length(a)] <- log(cutoff[[j]][3 : (length(cutoff[[j]]) - 1)] - cutoff[[j]][2 : (length(cutoff[[j]]) - 2)])
        aprop <- a[2 : length(a)] + rnorm(length(2 : length(a)), 0, sqrt(0.1))
        cutoffprop <- cutoff[[j]]
        cutoffprop[3] <- exp(aprop)[1]
        if (ncat[j] >= 4)
        {
          for (jord in 4 : ncat[j])
          {
            cutoffprop[jord] <- sum(exp(aprop)[1 : (jord - 2)])
          }
        }
        TUprop <- cutoffprop[Y[, j] + 2] - XBgamma
        TU[, j] <- cutoff[[j]][Y[, j] + 2] - XBgamma
        TLprop <- cutoffprop[Y[, j] + 1] - XBgamma
        TL[, j] <- cutoff[[j]][Y[, j] + 1] - XBgamma
        priorRatio <- sum(dnorm(aprop, 0, sqrt(10), log = TRUE)) - sum(dnorm(a[2 : length(a)], 0, sqrt(10), log = TRUE))
        DiffLikProp <- pnorm((TUprop[indna] - conmu[indna, j]) / sqrt(convar[indna, j])) - pnorm((TLprop[indna] - conmu[indna, j]) / sqrt(convar[indna, j]))
        DiffLikProp[which(DiffLikProp == 0)] < pnorm(Inf) - pnorm(8.2)
        
        ########## comment ##########
        # In equation (12) of the main paper (Alexopoulos and Bottolo, 2019), there is the difference of two standard Gaussian cdfs, 
        # but this difference can be come very small (zero in R) if the cdfs are both evaluated at large values. Instead, we compute 
        # pnorm(Inf) - pnorm(8.2) which is the smallest difference that can be computed in R
        
        logLikProp <- sum(log(DiffLikProp))
        DiffLik <- ((pnorm((TU[indna, j] - conmu[indna, j]) / sqrt(convar[indna, j])) - pnorm((TL[indna, j] - conmu[indna, j]) / sqrt(convar[indna, j]))))
        DiffLik[which(DiffLik == 0)] <- pnorm(Inf) - pnorm(8.2)
        
        ########## comment ##########
        # In equation (12) of the main paper (Alexopoulos and Bottolo, 2019), there is the difference of two standard Gaussian cdfs, 
        # but this difference can be come very small (zero in R) if the cdfs are both evaluated at large values. Instead, we compute 
        # pnorm(Inf) - pnorm(8.2) which is the smallest difference that can be computed in R
        
        logLikOld <- sum(log(DiffLik))
        likRatio <- logLikProp - logLikOld
        if (log(runif(1)) < (likRatio + priorRatio))
        {
          cutoff[[j]][3] <- cutoffprop[3]
          if (ncat[j] >= 4)
          {
            for (jord in 4 : ncat[j])
            {
              cutoff[[j]][jord] <- cutoffprop[jord]
            }
          }
          TU[, j] <- TUprop
          TL[, j] <- TLprop
        }
      }
      
      ########## update negative binomial parameter ##########
      
      if (response_type[j] == "negative binomial")
      {
        mygamma <- which(gammacurr[start : end] > 0)
        if (length(mygamma) > 0)
        {
          Xgamma <- as.matrix(Xactual[, mygamma])
          Bj <-Bcurr[start : end]
          Bjgamma <- Bj[mygamma]
          XBgamma <- Xgamma %*% Bjgamma
          transXgamma <-t(Xgamma)
        }else{
          Xgamma <- matrix(0, n, 1)
          Bj <- rep(0, commoncovs)
          Bjgamma <- 0
          XBgamma <- rep(0, n)
          transXgamma <- t(Xgamma)
        }
        formu <- Rcurr[j, -j] %*% solve(Rcurr[-j, -j])
        conmu[, j] <- formu %*% t(ZZ[, -j])
        convar[,j] <- Rcurr[j, j] - formu %*% Rcurr[-j, j]
        rhoprop <- log(negbinompar[myj]) + rnorm(1, 0, sqrt(0.05))
        TUprop <- qnorm(pnbinom(Y[, j], exp(rhoprop), 1 - (1 /(1 + exp(-XBgamma)))))
        TLprop <- qnorm(pnbinom(Y[, j] - 1, exp(rhoprop), 1 - (1 /(1 + exp(-XBgamma)))))
        yy <- union(which(TUprop == TLprop), which(TLprop == Inf))
        if ((length(yy) > 0))
        {
          TLprop[yy] <- qnorm(0.9999999999999)
          TUprop[yy] <- qnorm(0.99999999999999)
          
          ########## comment ##########
          # For some outliers $y_{ik}$, the pdf $F_{k}(y_{ik}|\beta_{k},\theta_{k})$ and $F_{k(y_{ik} - 1|\beta_{k},\theta_{k})$ can
          # become one in R and qnorm(1) = Inf, but TL and TU cannot be equal since they are the truncation points for the Gaussian 
          # distribution from which the latent variables are drawn. Therefore, whenever $F_{k}(y_{ik}|\beta_{k},\theta_{k})$=1 in
          # R, we set $F_{k}(y_{ik}|\beta_{k},\theta_{k})$=0.99999999999999 and $F_{k(y_{ik} - 1|\beta_{k},\theta_{k})$= 0.9999999999999
        
        }
        priorRatio <- dgamma(exp(rhoprop), 2, rate = 1, log = TRUE) - dgamma(negbinompar[myj], 2, rate = 1, log = TRUE)
        DiffLikProp <- pnorm((TUprop[indna] - conmu[indna, j]) / sqrt(convar[indna, j])) - pnorm((TLprop[indna] - conmu[indna, j]) / sqrt(convar[indna, j]))
        DiffLikProp[which(DiffLikProp == 0)] <- pnorm(Inf) - pnorm(8.2)
        
        ########## comment ##########
        # In equation (12) of the main paper (Alexopoulos and Bottolo, 2019), there is the difference of two standard Gaussian cdfs, 
        # but this difference can be come very small (zero in R) if the cdfs are both evaluated at large values. Instead, we compute 
        # pnorm(Inf) - pnorm(8.2) which is the smallest difference that can be computed in R
        
        logLikProp <- sum(log(DiffLikProp))
        DiffLik <- ((pnorm((TU[indna, j] - conmu[indna, j]) / sqrt(convar[indna, j])) - pnorm((TL[indna, j] - conmu[indna, j]) / sqrt(convar[indna, j]))))
        DiffLik[which(DiffLik == 0)] <- pnorm(Inf) - pnorm(8.2)
        
        ########## comment ##########
        # In equation (12) of the main paper (Alexopoulos and Bottolo, 2019), there is the difference of two standard Gaussian cdfs, 
        # but this difference can be come very small (zero in R) if the cdfs are both evaluated at large values. Instead, we compute 
        # pnorm(Inf) - pnorm(8.2) which is the smallest difference that can be computed in R
        
        logLikOld <- sum(log(DiffLik))
        likRatio <- logLikProp-logLikOld
        if (log(runif(1)) < (likRatio + priorRatio + rhoprop - log(negbinompar[myj])))
        {
          negbinompar[myj] <- exp(rhoprop)
          TU[, j] <- TUprop
          TL[, j] <- TLprop
        }
      }
      
      ########## update Z variables ##########
      
      for (i in 1 : n)
      {
        if (!is.na(Y[i, j]))
        {
          eps <- rtruncnorm(1, a = TL[i, j], b = TU[i, j], formu %*% (ZZ[i, -j]), sqrt(Rcurr[j, j] - formu %*% Rcurr[-j, j]))
          ZZ[i, j] <- eps
          if (is.na(ZZ[i, j]))
          {
            cat(c(TL[i, j], TU[i, j]), "\n")
          }
        }else{
          eps <- rnorm(1, formu %*% (ZZ[i, -j]), sqrt(Rcurr[j, j] - formu %*% Rcurr[-j, j]))
          ZZ[i, j] <- eps
        }
      }
    }
  }
  output <- list(gammacurr, Bcurr, XB, ZZ, cutoff, negbinompar)
    
  return(output)
}

#' @title BVS4GCR_MCMC
#' @description MCMC implementation of Bayesian Variable Selection for Gaussian Copula Regression models
#'
#' @param nits Number of MCMC iterations (including burn-in)
#' @param burnin Number of iterations to be discarded as burn-in
#' @param thin Store the outcome at every thin-th iteration
#' @param monitor Display the monitored-th iteration
#' @param Y (number of samples) times (number of responses) matrix of responses
#' @param X (number of samples) times (number of predictors) matrix of predictors. If the model contains an intercept (for all 
#' responses), the first column is a vector of ones, followed by the covariates and the predictors on which Bayesian Variable 
#' Selection is performed
#' @param pfixed Number of covariates (including the intercept)
#' @param response_type Character vector specifying the response's type. Response's types currently supported are \code{c("Gaussian", 
#' "binary", "ordinal", "binomial", "negative binomial")}
#' @param beta_binom Parameters of the beta-binomial prior distribution. For details, see \insertCite{Kohn2001;textual}{BVS4GCR} and
#' \insertCite{Alexopoulos2019;textual}{BVS4GCR}
#' @param missing_cont Logical parameter that indicates the presence of missing values in any continuous response
#' @param fullcov Logical parameter to select full or sparse inverse correlation
#' @param nondecomp Logical, \code{TRUE} if a non-decomposable graphical model is selected
#' @param ntrials Number of trials of binomial responses if \code{"binomial"} is included in \code{response_type}
#' @param negbin_init Initial value for the over-dispersion parameter of the negative binomial distribution
#'
#' @details
#' For details regarding the model and the algorithm, see details in \insertCite{Alexopoulos2019;textual}{BVS4GCR}. Types of 
#' responses currently supported: \code{c("Gaussian", "binary", "ordinal", "binomial", "negative binomial")}. For ordinal categorical 
#' responses it is assumed that the first category is labelled with zero. The maximum number of categories supported is 5.
#'
#' @export
#'
#' @return The value returned is a list object \code{list(Gamma, B, G, R, D, Cutoff, NGB)}
#' \itemize{
#'   \item{\code{Β}}{ matrix of the (thinned) samples drawn from the posterior distribution of the regression coefficients}
#'   \item{\code{G}}{ 3D array of the (thinned) samples drawn from the posterior distribution of the graphs}
#'   \item{\code{R}}{ 3D array of the (thinned) samples drawn from the posterior distribution of the correlation matrix between the
#'   responses}
#'   \item{\code{D}}{ matrix of the (thinned) samples drawn from the posterior distribution of the standard deviations of the 
#'   responses (non-identifiable parameters in the case of discrete responses)}
#'   \item{\code{Cutoff}}{ list of the samples drawn from the posterior distribution of the cut-off points for ordinal categorical 
#'   responses}
#'   \item{\code{OD}}{ list of the samples drawn from the posterior distribution of the over-dispersion parameter for the negative 
#    binomial distributions} 
#'   responses}
#' }
#'
#' @references
#'   \insertAllCited{}
#'
#' @examples
### 100 characters ################################################################################
#' # Example 1: Toy example with only Gaussian responses and decomposable graphs specification
#'
#' set.seed(12345)
#' d <- Sim_GCRdata(4, 500, 100, 2, rep("Gaussian", 4), extras = list(rhoY = 0.8, rhoX = 0.7, 
#' fixed_effect = rep(1, 2), mu_effect = rep(3, 4), var_effect = rep(0.5, 4), s1Lev = 0.05, 
#' s2Lev = 0.95, std_responses = rep(1, 4)))
#'
#' output <- BVS4GCR_MCMC(250, 50, 4, 100, d$Y, d$X, 2, rep("Gaussian", 4), c(5, 3 ^2))
#'
#'
#' # Example 2: Toy example with only Gaussian responses and non-decomposable graphs 
#'
#' set.seed(12345)
#' d <- Sim_GCRdata(8, 1000, 100, 1, rep("Gaussian", 8), extras = list(rhoY = 0.8, rhoX = 0.7, 
#' fixed_effect = 1, mu_effect = rep(3, 8), var_effect = rep(0.5, 8), s1Lev = 0.05, s2Lev = 0.95, 
#' std_responses = rep(1, 8)))
#'
#' output <- BVS4GCR_MCMC(250, 50, 4, 100, d$Y, d$X, 1, rep("Gaussian", 8), c(5, 3 ^2), 
#' nondecomp = TRUE)
#'
#'
#' # Example 3: Combination of Gaussian, binary and ordinal responses (similar to Scenario I in
#' # Alexopoulos and Bottolo (2019))
#'
#' set.seed(12345)
#' d <- Sim_GCRdata(4, 50, 30, 1, c("Gaussian", "binary", "ordinal", "ordinal"), extras = 
#' list(rhoY = 0.8, rhoX = 0.7, fixed_effect = 1, mu_effect = c(1, rep(0.5, 3)), var_effect = 
#' c(1, rep(0.2, 3)), s1Lev = 0.15, s2Lev = 0.95, std_responses = rep(1, 4), ncat = c(3, 4), 
#' cutpoints = cbind(c(0.5, NA, NA), c(0.5, 1.2, 2))))
#'
#' output <- BVS4GCR_MCMC(250, 50, 4, 100, d$Y, d$X, 1, c("Gaussian", "binary", "ordinal", 
#' "ordinal"), c(5, 3 ^2))
#'
#'
#' # Example 4: Combination of Gaussian and count responses (Scenario IV in Alexopoulos and Bottolo 
#' # (2019))
#'
#' set.seed(12345)
#' d <- Sim_GCRdata(4, 100, 100, 2, c("Gaussian", "binomial", "negative binomial", 
#' "negative binomial"), extras = list(rhoY = 0.8, rhoX = 0.7, fixed_effect = c(-0.5, 0.5), 
#' mu_effect = c(1, rep(0.5, 3)), var_effect = c(1, rep(0.2, 3)), s1Lev = 0.05, s2Lev = 0.95, 
#' std_responses = rep(1, 4), binom_par = 10, negbin_par = c(0.5, 0.75)))
#'
#' output <- BVS4GCR_MCMC(250, 50, 4, 100, d$Y, d$X, 2, c("Gaussian", "binomial", 
#' "negative binomial", "negative binomial"), c(5, 3 ^2), ntrials = 10, negbin_init = 
#' c(0.5, 0.75))
#'
#' @importFrom stats binomial
#' @importFrom stats coefficients
#' @importFrom stats cor
#' @importFrom stats cov2cor
#' @importFrom stats dbinom
#' @importFrom stats dgamma
#' @importFrom stats dgeom
#' @importFrom stats dnorm
#' @importFrom stats glm
#' @importFrom stats lm
#' @importFrom stats pbinom
#' @importFrom stats pgeom
#' @importFrom stats pnbinom
#' @importFrom stats pnorm
#' @importFrom stats qgeom
#' @importFrom stats qnorm
#' @importFrom stats rbinom
#' @importFrom stats rchisq
#' @importFrom stats rgamma
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats var
#' @importFrom MASS mvrnorm
#' @importFrom MASS polr
#' @importFrom MCMCpack rwish
#' @importFrom mvnfast rmvn
#' @importFrom mvtnorm dmvnorm
#' @importFrom truncnorm rtruncnorm
#' @importFrom Rdpack reprompt


BVS4GCR_MCMC <- function(nits, burnin, thin, monitor, Y, X, pfixed, response_type, beta_binom = c(3, 1), 
                         missing_cont = FALSE, fullcov = FALSE, nondecomp = FALSE, 
                         ntrials = NULL, negbin_init = NULL)
{
  # print(nits)
  # print(burnin)
  # print(thin)
  # print(monitor)
  # print(pfixed)
  # print(response_type)
  # print(beta_binom)
  # print(missing_cont)
  # print(missing_cont)
  # print(nondecomp)
  # print(ntrials)
  # print(negbin_init)
  
  n <- dim(Y)[1]
  m <- dim(Y)[2]
  p <- dim(X)[2]
  surecovs <- pfixed
  ncat <- rep(0, m)
  negbinompar <- negbin_init
  
  ########## compute pvalues needed for the proposal for gamma ##########
  
  pvalues <- cors <- matrix(NA, p - pfixed, m)
  ii <- 0
  for (i in 1 : (p - pfixed))
  
  ########## comment ##########
  # In the next code release, "pfixed" covariates will be included in the GLM model to inform pvalues more correctly
  
  {
    for (j in 1 : m)
    {
      myj <- 0
      myj2 <- 0
      cors[i, j] <- abs(cor(Y[, j], X[, i + surecovs], use = "complete.obs"))
      if (response_type[j] == "Gaussian")
      {
        pvalues[i, j] <- summary(lm(Y[, j] ~ X[, i + surecovs]))[[4]][2, 1]
      }else{
        thelink = "probit"
        
        ########## comment ##########
        # Hard coding! probit link is used
        
        if (response_type[j] == "binary")
        {
          fit <- glm(Y[, j] ~ X[, i + surecovs], family = binomial(link = thelink))
          pvalues[i, j] <- summary(fit)[[12]][2, 1]
        }
        if (response_type[j] == "ordinal")
        {
          pvalues[i,j] <- coefficients(polr(as.factor(Y[,j]) ~ X[, i + surecovs], method = thelink))[[1]]
        }
        if (response_type[j] == "binomial")
        {
          myj <- myj + 1
          fit <- glm((Y[, j] / ntrials[myj]) ~ X[, i + surecovs], family = binomial(link = "logit"), weights = rep(ntrials[myj], n))
                 pvalues[i, j] <- summary(fit)[[12]][2, 4]
        }
        if (response_type[j] == "negative binomial")
        {
          fit <- glm(Y[, j] ~ X[, i + surecovs], family = "poisson")
                 pvalues[i, j] <- summary(fit)[[12]][2,4]
        }
      }
    }
  }
  cors <- abs(pvalues)
  
  ########## beta-binomial prior on gamma ##########
  
  Eqg <- beta_binom[1]
  Varqg <- beta_binom[2]
  K <- p - pfixed
  Epi <- Eqg / K
  meanbeta <- Epi
  varbeta <- (Varqg - K * Epi * (1 - K * Epi)) / (K * (K - 1) * Epi)
  api <- (meanbeta * varbeta - meanbeta) / (meanbeta - varbeta)
  bpi <- api * (1 - meanbeta) / meanbeta
  
  mc <- mo <- mordinal <- mbinary <- mbinom <- mnegbinom <- 0
  for (i in 1 : length(response_type))
  {
    if (response_type[i] == "Gaussian")
    {
      mc <- mc + 1
    }
    if (response_type[i] == "binary")
    {
      mbinary <- mbinary + 1
      mo <- mo + 1
    }
    if (response_type[i] == "ordinal")
    {
      mordinal <- mordinal + 1
      mo <- mo+1
    }
    if (response_type[i] == "binomial")
    {
      mbinom <- mbinom + 1
      mo <- mo+1
    }
    if (response_type[i] == "negative binomial")
    {
      mnegbinom = mnegbinom + 1
      mo = mo + 1
    }
  }
  
  ########## starting values for gamma ##########
  
  gammacurr <- rep(0, p * m)
  inds <- rep(NA, m * surecovs)
  if (surecovs > 0)
  {
    for(kk in c(0, 1 : (m - 1)))
    {
      inds[(kk * surecovs + 1) : ((kk + 1) * surecovs)] <- (1 : surecovs) + (kk * (p - pfixed) + kk * surecovs)
    }
    gammacurr[inds] <- 1
  }
  for (kk in c(0, 1 : (m - 1)))
  {
    thres <- mean(cors[, kk + 1]) + 3 * sqrt(var(cors[, kk + 1]))
    gammacurr[which(cors[, kk + 1] >= thres) + (kk * (p - pfixed) + (kk + 1) * surecovs)] <- 1
  }
  
  ########## create matrices to store MCMC output ##########
  
  Bsave <- matrix(NA, (nits - burnin) /thin, p * m)
  Gsave <- Rsave <- array(NA, c(m, m, (nits - burnin) / thin))
  Dsave <- matrix(NA, m, (nits - burnin) / thin)
  cutoffsave <- list()
  negbinomparsave <- list()
  jj = 0
  
  ########## initial values for decomposable graph, correlation matrix and betas ##########
  
  Gcurr <- diag(rep(1, m))
  diag(Gcurr) <- 0
  if ((mbinary > 0) || (mordinal > 0))
  {
    cutoffproplow <- cutoffpropup <- matrix(NA, n, mbinary + mordinal)
  }
  Rcurr <- diag(rep(1, m))
  invRcurr <- solve(Rcurr)
  Bcurr <- rep(0, p * m)
  Bcurr[which(gammacurr != 0)] <- rnorm(length(which(gammacurr != 0)), 0, 0.1)
  myB <- matrix(Bcurr, m, p, byrow = TRUE)
  d <- rep(1, m)
  D <- diag((d))
  XB <- matrix(NA, n, m)
  for (i in 1 : n)
  {
    for(k in 1 : m)
    {
      XB[i, k] = sum(X[i,] * myB[k, ])
    }
  }
  invSigma <- diag(1 /d) %*% invRcurr %*% diag(1 /d)
  Sigma <- solve(invSigma)
  
  ########## initial values for Z ##########
  
  Zcurr <- matrix(0, n, m)
  if (mc > 0)
  {
    Zcurr[, which(response_type == "Gaussian")] <- (Y[, which(response_type == "Gaussian")] - XB[, which(response_type == "Gaussian")])
    Zcurr[, which(response_type == "Gaussian")] <- Zcurr[, which(response_type == "Gaussian")] / sqrt(d[which(response_type == "Gaussian")])
  }
  if (missing_cont == TRUE)
  {
    Zcurr[, 1 : mc][which(is.na(Zcurr[, 1 : mc]))] <- rnorm(1)
  }
  
  ########## initial values for cut-point and delta ##########
  
  pi_delta <- 10 ^(-4)
  v1 <- 1
  v0 <- 0.01
  deltas <- matrix(0, m, m)
  vv <- matrix(0.01, m, m)
  cutoff <- list()
  for (k in 1 : m)
  {
    if ((response_type[k] == "Gaussian") | (response_type[k] == "binomial") | (response_type[k] == "negative binomial"))
    {
      cutoff[[k]] <- NA
    }
    if (response_type[k] == "binary")
    {
      cutoff[[k]] <- c(-Inf, 0, Inf)
    }
    if (response_type[k] == "ordinal")
    {
      nolvs <- max(Y[, k]) + 1
      ncat[k] <- nolvs
      cutoff[[k]] <- rep(NA, nolvs + 1)
      cutoff[[k]][c(1, 2, nolvs + 1)] <- c(-Inf, 0, Inf)
      jjjj <- 0
      for(jjj in 3 : nolvs)
      {
        cutoff[[k]][jjj] <- 0.5 + jjjj
        jjjj <- jjjj + 1.5
      }
    }
  }
  Wcurr <- Zcurr %*% D
  
  ########## initial values for non-decomposable graph ##########
  
  if (nondecomp)
  {
    gg <- mod_BDgraph_init(Wcurr, niter = 100)
    
    ########## comment ##########
    # Hard coding! "iter" is the internal number of iterations setup to generate an initial non-decomposable graph. Only the last 
    # simulated graphical model is retained 
    
    invSigma <- gg$last_K
    Sigma <- solve(invSigma)
    d <- sqrt(diag(Sigma))
    Rcurr <- diag(1 /d) %*% Sigma %*% diag(1 /d)
    invRcurr <- diag(d) %*% invSigma %*% diag(d)
    D <- diag((d))
    Zcurr <- Wcurr %*% diag(1/d)
  }else{
    Gcurr = diag(m)
  }
  
  ########## start MCMC ########## 
  
  for (j in 1 : nits)
  {
    Latent <- BVS(cors, Y, Bcurr, XB, X, Rcurr, invRcurr, Sigma, invSigma, Zcurr, cutoff, 
                  response_type, api, bpi, ncat, p - pfixed, gammacurr, surecovs, d, negbinompar, ntrials)
    
    gammacurr <- Latent[[1]]
    Bcurr <- Latent[[2]]
    XB <- Latent[[3]]
    Zcurr <- Latent[[4]]
    cutoff <- Latent[[5]]
    negbinompar <- Latent[[6]]
    Wcurr <- Zcurr %*% D
    
    if (fullcov == FALSE)
    {
      if (nondecomp)
      {
        gg <- mod_BDgraph(data = Wcurr, niter = 100, burnin = 30, df.prior = 2, g.start = gg$last_graph, K.start = gg$last_K, cores = 1)
        
        ########## comment ##########
        # Hard coding! "niter" is the internal number of iterations setup to simulate non-decomposable graphs, "burnin" is the number 
        # of burn-in iterations (overall iterations are niter - burnin) and df.prior are the prior degree of freedom. Only the last 
        # simulated graphical model is retained 
        
        invSigma <- gg$last_K
        Sigma <- solve(invSigma + diag(rep(10 ^(-10)), m))
        Gcurr <- gg$last_graph
      }else{
        constrG <- Sample_G(Wcurr, Gcurr)
        Gcurr <- constrG
        rhiwish <- Sim_HIW(Gcurr, t(Wcurr) %*% Wcurr + diag(rep(1, m)), n + 2)
        Sigma <- rhiwish[[1]]
        invSigma <- rhiwish[[2]]
      }
    }
    if (fullcov == TRUE)
    {
      invSigma <- rwish(1 + n + m, solve(t(Wcurr) %*% Wcurr + diag(rep(1 ,m))))
      Sigma <- solve(invSigma + diag(rep(10 ^(-10), m)))
    }
    d <- sqrt(diag(Sigma))
    Rcurr <- diag(1 /d) %*% Sigma %*% diag(1 /d)
    invRcurr <- diag(d) %*% invSigma %*% diag(d)
    D <- diag((d))
    Zcurr <- Wcurr %*% diag(1 /d)
    if (j > burnin)
    {
      if (j %% thin == 0)
      {
        jj <- jj + 1
        Bsave[jj, ] <- Bcurr
        Gsave[, , jj] <- Gcurr
        Rsave[, , jj] <- Rcurr
        Dsave[, jj] <- d
        
        cutoffsave[[jj]] <- cutoff
        negbinomparsave[[jj]] <- negbinompar
      }
    }
    if (j %% monitor == 0)
    {
      cat(paste("iteration", j), "\n")
    }
  }
  output <- list(B = Bsave, G = Gsave, R = Rsave, D = Dsave, Cutoff = cutoffsave, OD = negbinomparsave)
  
  return(output)
}

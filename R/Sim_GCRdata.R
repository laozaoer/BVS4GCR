#' @title Sim_GCRdata
#'
#' @description Simulate data from the Gaussian Copula Regression model
#'
#' @param m Integer, number of responses
#' @param n Integer, number of samples
#' @param p Integer, number of predictors
#' @param pfixed Integer, number of covariates (including the intercept). If \code{pfixed = 1} only the intercept will be simulated
#' @param response_type \code{m}-dimensional character vector specifying the response's type. Response's types currently supported are 
#' \code{c("Gaussian", "binary", "ordinal", "binomial", "negative binomial")}
#' @param extras List specifying response-specific arguments. These are the one-lag correlation \code{rhoY} in the autoregressive 
#' model used to simulate the responses, the correlation \code{rhoX} between the predictors (for both, see details in 
#' \insertCite{Rothman2010;textual}{BVS4GCR}), \code{fixed_effect} the vector of the regression coefficients for the fixed 
#' predictors (if simulated) where the first element is assumed to be the intercept, the mean \code{mu_effect} and variance \code{var_effect} 
#' of the Gaussian distribution used to simulate the regression coefficients, the probabilities \code{s1Lev} and \code{s2Lev} utilized 
#' to induce sparsity on the regression coefficients, the standard deviation \code{std_responses} of the response variables, the number
#' of categories \code{ncat} of each ordinal categorical response (if simulated), the cut-off points \code{cutpoints} used to simulate 
#' the ordinal variables (if simulated), \code{binom_par} and \code{negbinom} the numer of trials in the binomial experiment and the 
#' over-dispersion parameter of the negative binomial distribution. See also the examples below and details in \insertCite{Alexopoulos2019;textual}{BVS4GCR}
#'
#' @details
#' The parameters \code{s1Lev} and \code{s2Lev} in the \code{extras} list can be specified following the fact that ((1 - \code{s2Lev}) * \code{p}) 
#' predictors are irrelevant for all \code{m} responses and that each relevant predictor is common across (\code{s1Lev} * \code{m}) 
#' predictors \insertCite{Rothman2010}{BVS4GCR}. Type of responses currently supported: (continuous) Gaussian, (discrete) binary,
#' ordinal categorical and count (binomial and negative binomial distribution)
#'
#' @export
#'
#' @return The value returned is a list object \code{list(Y, X, B, R)}
#' \itemize{
#'  \item{\code{Y}}{ \code{n} times \code{m} matrix of responses}
#'  \item{\code{X}}{ \code{n} times (\code{pfixed} + \code{p}) matrix of covariates (including the intercept if simulated) and 
#'  predictors on which Bayesian Variable 
#'  Selection is performed}
#'  \item{\code{B}}{ \code{m} times (\code{pfixed} + \code{p}) matrix of regression coefficients}
#'  \item{\code{R}}{ \code{m} times \code{m} correlation matrix used to simulate the responses}
#' }
#'
#' @references
#' \insertAllCited{}
#'
### 100 characters ################################################################################
#' @examples
#' # Example 1: Simulate a combination of Gaussian, binary and ordinal responses
#'
#' set.seed(12345)
#' d <- Sim_GCRdata(4, 500, 30, 1, c("Gaussian", "binary", "ordinal", "ordinal"), extras = 
#' list(rhoY = 0.8, rhoX = 0.7, fixed_effect = 1, mu_effect = c(1, rep(0.5, 3)), var_effect = 
#' c(1, rep(0.2, 3)), s1Lev = 0.15, s2Lev = 0.95, std_responses = rep(1, 4), ncat = c(3, 4), 
#' cutpoints = cbind(c(0.5, NA, NA), c(0.5, 1.2, 2))))
#'
#' par(mfrow = c(2, 2))
#' hist(d$Y[,1], main = "Gaussian", xlab = expression(Y[1]), prob = TRUE, ylab = "")
#' barplot(prop.table(table(d$Y[, 2])), main = "binary", xlab = expression(Y[2]))
#' barplot(prop.table(table(d$Y[, 3])), main = "ordinal (3 categories)", xlab = expression(Y[3]))
#' barplot(prop.table(table(d$Y[, 4])), main = "ordinal (4 categories)", xlab = expression(Y[4]))
#' 
#' par(mfrow = c(2, 2))
#' plot(d$B[1, ], ylab=expression(beta[1]), main=expression(Y[1]), xlab="Pred.",pch = 16,cex = .7)
#' plot(d$B[2, ], ylab=expression(beta[2]), main=expression(Y[2]), xlab="Pred.",pch = 16,cex = .7)
#' plot(d$B[3, ], ylab=expression(beta[3]), main=expression(Y[3]), xlab="Pred.",pch = 16,cex = .7)
#' plot(d$B[4, ], ylab=expression(beta[4]), main=expression(Y[4]), xlab="Pred.",pch = 16,cex = .7)
#'
#'
#' # Example 2: As in Example 1 with more categories for the forth response 
#'
#' set.seed(12345)
#' d <- Sim_GCRdata(4, 500, 30, 1, c("Gaussian", "binary", "ordinal", "ordinal"), extras = 
#' list(rhoY = 0.8, rhoX = 0.7, fixed_effect = 1, mu_effect = c(1, rep(0.5, 3)), var_effect = 
#' c(1, rep(0.2, 3)), s1Lev = 0.15, s2Lev = 0.95, std_responses = rep(1, 4), ncat = c(3, 5),
#' cutpoints = cbind(c(0.5, 1, NA, NA), c(0.5, 0.7, 1.5, 2))))
#'
#' par(mfrow = c(2, 2))
#' hist(d$Y[, 1], main = "Gaussian", xlab = expression(Y[1]), prob = TRUE, ylab = "")
#' barplot(prop.table(table(d$Y[, 2])), main = "binary", xlab = expression(Y[2]))
#' barplot(prop.table(table(d$Y[, 3])), main = "ordinal (3 categories)", xlab = expression(Y[3]))
#' barplot(prop.table(table(d$Y[, 4])), main = "ordinal (5 categories)", xlab = expression(Y[4]))
#'
#'
#' # Example 3: Simulate a combination of Gaussian and count responses
#'
#' set.seed(12345)
#' d <- Sim_GCRdata(4, 1000, 100, 2, c("Gaussian", "binomial", "negative binomial", 
#' "negative binomial"), extras = list(rhoY = 0.8, rhoX = 0.7, fixed_effect = c(-0.5, 0.5), 
#' mu_effect = c(1, rep(0.5, 3)), var_effect = c(1, rep(0.2, 3)), s1Lev = 0.05, s2Lev = 0.95, 
#' std_responses = rep(1, 4), binom_par = 10, negbin_par = c(0.5, 0.75)))
#'
#' par(mfrow = c(2, 2))
#' hist(d$Y[,1], main = "Gaussian", xlab = expression(Y[1]), prob = TRUE, ylab = "")
#' hist(table(d$Y[, 2]), main = "binomial", xlab = expression(Y[2]))
#' hist(table(d$Y[, 3]), main = "negative binomial", xlab = expression(Y[3]))
#' hist(table(d$Y[, 4]), main = "negative binomial", xlab = expression(Y[4]))
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


Sim_GCRdata <- function(m, n, p, pfixed, response_type, 
                        extras = list(rhoY = 0.8, rhoX = 0.7, fixed_effect = rep(1, pfixed), mu_effect = rep(0, m), var_effect = rep(1, m), 
                        s1Lev = 0.10, s2Lev = 1, std_responses = rep(1, m), 
                        ncat = NULL, cutpoints = NULL, binom_par = NULL, negbin_par = NULL))
{
  # print(m)
  # print(n)
  # print(p)
  # print(pfixed)
  # print(response_type)
  # print(extras)
  
  p <- p + pfixed
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
      mo <- mo + 1
    }
    if (response_type[i] == "binomial")
    {
      mbinom <- mbinom + 1
      mo <- mo + 1
    }
    if (response_type[i] == "negative binomial")
    {
      mnegbinom <- mnegbinom + 1
      mo <- mo + 1
    }
  }
  Z <- simY <- matrix(NA, n, m)
  
  ########## create correlation matrix for responses ##########
  
  Sigma <- matrix(NA, m, m)
  for (i in 1 : m)
  {
    for (j in 1 : m)
    {
      Sigma[i, j] <- extras$rhoY ^(abs(i - j))
    }
  }
  
  ########## create band covariance for predictors ##########
  
  corsX <- matrix(0, p - pfixed, p - pfixed)
  for (i in 1 : (p - pfixed))
  {
    for (j in 1 : (p - pfixed))
    {
      corsX[i, j] <- extras$rhoX ^(abs(i - j))
    }
  }
  
  ########## simulate predictors ##########
  
  simX <- matrix(NA, n, p)
  if (pfixed >= 1)
  {
    simX[, 1] <- rep(1, n)
    if (pfixed >= 2)
    {
      simX[, 2 : pfixed] <- rnorm(n * length(2 : pfixed))
    }
  }
  for (i in 1 : n)
  {
    simX[i, (pfixed + 1) : p] <- rmvn(1, rep(0, p - pfixed), corsX)
  }
  
  ########## simulate effects ##########
  
  commoncovariates <- p - pfixed
  K <- matrix(0, commoncovariates, m)
  W <- matrix(0, commoncovariates, m)
  for (k in 1 : m)
  {
    K[, k] <- rbinom(commoncovariates, 1, extras$s1Lev)
    W[, k] <- rnorm(commoncovariates, extras$mu_effect[k], sqrt(extras$var_effect[k]))
  }
  Q <- matrix(0, commoncovariates, m)
  rows_ind <-rbinom(commoncovariates, 1, extras$s2Lev)
  Q[which(rows_ind > 0), ] <-1
  BB <- t(W * K * Q)
  simB <- matrix(0, m, p)
  if (pfixed > 0)
  {
    fixed_effect <- extras$fixed_effect
    for (k in 1 : m)
    {
      simB[k, c(1 : pfixed, which(BB[k,] != 0) + pfixed)] <- c(fixed_effect, BB[k, c(which(BB[k, ] != 0))])
    }
  }else{
    for (k in 1 : m)
    {
      simB[k, c(which(BB[k, ] != 0))] <- c(BB[k, c(which(BB[k, ] != 0))])
    }
  }
  Bvec <- as.vector(t(simB))
  R <- cov2cor(Sigma)
  D <- diag(extras$std_responses)
  R <- D %*% R %*% D
  XB <- matrix(NA, n, m)
  for (i in 1 : n)
  {
    Z[i, ] <- mvrnorm(1, rep(0, m), R)
    for (k in 1 : m)
    {
      XB[i,k] <- sum(simX[i, ] * simB[k, ])
    }
  }
  
  ########## simulate Gaussian responses ##########
  
  if (mc > 0)
  {
    simY[, 1 : mc] <- Z[, 1 : mc] + XB[, 1 : mc]
  }
  
  ########## simulate binary, categorical, count responses by specifying suitable cdfs ##########
  
  ########## comment ##########
  # The maximum number of categories is limited to 5. In the next code release, "ncat" will be extended to any number of categories.
  
  if (mo > 0)
  {
    if (mbinary > 0)
    {
      for (i in 1 : n)
      {
        eps <- rnorm(1, 0, sqrt(1))
        for (j in (mc + 1) : (mc + mbinary))
        {
          if ((Z[i, j] > -Inf) & (Z[i, j] <= -XB[i, j]))
          {
            simY[i, j] <- 0
          }
          if ((Z[i, j] > -XB[i, j]) & (Z[i, j] <= Inf))
          {
            simY[i, j] <- 1
          }
        }
      }
    }
    
    ########## ordinal categorical responses ##########
    
    if (mordinal > 0)
    {
      jj <- 0
      for (j in (mc + mbinary + 1) : (mc + mbinary + mordinal))
      {
        jj <- jj + 1
        ncat <- extras$ncat[jj]
        for (i in 1 : n)
        {
          if ((Z[i, j] <= -XB[i, j]) & (Z[i, j] >- Inf))
          {
            simY[i, j] <- 0
          }
          if ((Z[i, j] > -XB[i, j]) & (Z[i, j] <= extras$cutpoints[1, jj] - XB[i, j]))
          {
            simY[i, j] <- 1
          }
          if (ncat == 3)
          {
            if ((Z[i, j] > extras$cutpoints[1, jj] - XB[i, j]))
            {
              simY[i, j] <- 2
            }
          }else{
            if ((Z[i, j] > extras$cutpoints[1, jj] - XB[i, j]) & (Z[i, j] <= extras$cutpoints[2, jj] - XB[i, j]))
            {
              simY[i, j] <- 2
            }
            if (ncat == 4)
            {
              if ((Z[i, j] > extras$cutpoints[2, jj] - XB[i, j]))
              {
                simY[i, j] <- 3
              }
            }else{
              if ((Z[i, j] > extras$cutpoints[2, jj] - XB[i, j] ) & (Z[i, j] <= extras$cutpoints[3, jj] - XB[i,j]))
              {
                simY[i, j] <- 3
              }
              if (ncat == 5)
              {
                if ((Z[i, j] > extras$cutpoints[3, jj] - XB[i, j]))
                {
                  simY[i, j] <- 4
                }
              }
            }
          }
        }
      }
    }
    
    ########## binomial responses ##########
    
    if (mbinom > 0)
    {
      simBinom_par <- extras$binom_par
      jbinom <- 0
      for (j in (mc + mbinary + mordinal + 1) : (mc + mbinary + mordinal + mbinom))
      {
        jbinom <- jbinom + 1
        for (i in 1: n)
        {
          k <- 0
          while((Z[i, j] > qnorm(pbinom(k, simBinom_par[jbinom], exp(XB[i, j]) / (1 + exp(XB[i,j]))))) | 
                (Z[i, j] < qnorm(pbinom(k - 1, simBinom_par[jbinom], exp(XB[i, j]) / (1 + exp(XB[i, j]))))))
          {
            k <- k+1
            # cat(k, "\n")
          }
          simY[i, j] <- k
        }
      }
    }
    
    ########## negative binomial responses ##########
    
    if (mnegbinom > 0)
    {
      negbin_par <- extras$negbin_par
      jnegbinom <- 0
      for (j in (mc + mbinary + mordinal + mbinom + 1) : (mc + mbinary + mordinal + mbinom + mnegbinom))
      {
        jnegbinom <- jnegbinom+1
        for (i in 1 : n)
        {
          k <- 0
          while ((Z[i, j] > qnorm(pnbinom(k, negbin_par[jnegbinom], 1 - (exp(XB[i, j]) / (1 + exp(XB[i, j])))))) | 
                 (Z[i, j] < qnorm(pnbinom(k - 1, negbin_par[jnegbinom], 1 - (exp(XB[i, j]) / (1 + exp(XB[i, j])))))))
          {
            k <- k + 1
            # cat(k, "\n")
          }
          simY[i, j] <- k
        }
      }
    }
  }
  output <- list(Y = simY, X = simX, B = simB, R = Sigma)
  
  return(output)
}

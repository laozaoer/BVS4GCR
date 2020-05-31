
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BVS4GCR: Bayesian Variable Selection for Gaussian Copula Regression models

This package performs Bayesian predictors and covariance selection for combinations of continuous and/or discrete responses using a recently proposed Bayesian Variable Selection (BVS) algorithm for Gaussian Copula Regression (GCR) models 
([Alexopoulos and Bottolo, 2019](https://arxiv.org/pdf/1907.08245.pdf)).

The current version of the **BVS4GCR** R-package allows the analysis of (continuous) Gaussian and (discrete) binary, ordinal categorical and count responses (modelled with a binomial or negative binomial distribution) as well as the specification of a decomposable or non-decomposable graphical model for the conditional dependence structure between the responses. The latter instance of the algorithm takes full advantage of the C++ computational efficiency of the recently proposed R-package BDgraph ([Mohammadi and Wit, 2019](https://www.jstatsoft.org/article/view/v089i03)) to sample non-decomposable graphs. The **BVS4GCR** model also allows for covariates and confounders to be included which are not subject to variable selection and always included in the model. The current version of the package can also analyse data with missing at random in both continuous and discrete responses. 

The algorithm returns the posterior samples of the non-zero regression coefficients, the visited graphs as well as the residual covariance matrix between the responses and the response-specific parameters, i.e., the cut-points of ordinal categorical responses or the over-dispersion parameter for count data modelled with a negative binomial distribution. 

To favour the comparison with other (linear and non-linear)Bayesian Variable Selection algorithms for single-response regression, **BVS4GCR** R-package also includes a function to simulate the examples presented in [Alexopoulos and Bottolo, 2019](https://arxiv.org/pdf/1907.08245.pdf).

## Installation

**BVS4GCR** has several dependencies [**BDgraph**](https://cran.r-project.org/web/packages/BDgraph/index.html), [**MASS**](https://cran.r-project.org/web/packages/MASS/index.html), [**MCMCpack**](https://cran.r-project.org/web/packages/MCMCpack/), [**mvnfast**](https://cran.r-project.org/web/packages/mvnfast/index.html), [**mvtnorm**](https://cran.r-project.org/web/packages/mvtnorm/index.html) and [**truncnorm**](https://cran.r-project.org/web/packages/truncnorm/index.html) that are automatically imported and installed.

The installation of **BVS4GCR** requires the following steps:

1.  Install the [**devtools**](https://github.com/r-lib/devtools)
    package. This can be done from
    [**CRAN**](https://cran.r-project.org/). Invoke R and then type
    
    ``` 
     install.packages('devtools')
    ```

2.  Load the **devtools** package
    
    ``` 
     library('devtools')
    ```

3.  Install **BVS4GCR** package by
    typing
    
    ``` 
     devtools::install_github('lb664/BVS4GCR')
    ```

4.  Finally, load the **BVS4GCR** package
    
    ``` 
     library('BVS4GCR')
    ```

## Example 1

Installation can be checked by running the following toy example with Gaussian responses. Data is simulated by using the following command

    set.seed(12345) 
    d <- Sim_GCRdata(4, 500, 100, 2, rep("Gaussian", 4), extras = list(rhoY = 0.8, rhoX = 0.7, 
    fixed_effect = rep(0.5, 2), mu_effect = rep(3, 4), var_effect = rep(0.5, 4), s1Lev = 0.05, s2Lev = 0.95, std_responses = rep(1, 4)))


The function generates 500 observations for four Gaussian responses with 100 predictors and two covariates common to all responses (including the intercept). The underlying correlation between responses and collinearity amongst the predictors is controlled by the parameters `rhoY` and `rhoX`, respectively, and the level of sparsity by `s1Lev` and `s2Lev`. For details regarding these parameters, see [Rothman et al. (2010)](https://www.tandfonline.com/doi/abs/10.1198/jcgs.2010.09188).

**BVS4GCR** can be run on the simulated data as follows
    
    output <- BVS4GCR_MCMC(250, 50, 4, 100, d$Y, d$X, 2, rep("Gaussian", 4), c(5, 3 ^2))

where (250 - 50) posterior samples of the quantities of interest are recorded at every four MCMC iterations after burn-in, with the output of the algorithm monitored every 100 iterations. The _a priori_ sparsity prior `c(5, 3 ^2)` implies a range of associated predictors between 0 and (5 + 3 * 3) for each response.

## Example 2

This example is similar to the first simulation experiment presented in [Alexopoulos and Bottolo (2019)](https://arxiv.org/pdf/1907.08245.pdf) where 50 observations are simulated from a combination of one Gaussian, one binary and two ordinal categorical responses

    set.seed(12345)
    d <- Sim_GCRdata(4, 50, 30, 1, c("Gaussian", "binary", "ordinal", "ordinal"), extras = list(rhoY = 0.8, rhoX = 0.7, 
    fixed_effect = 1, mu_effect = c(1, rep(0.5, 3)), var_effect = c(1, rep(0.2, 3)), s1Lev = 0.15, s2Lev = 0.95, std_responses = rep(1, 4), ncat = c(3, 4), cutpoints = cbind(c(0.5, NA, NA), c(0.5, 1.2, 2))))

The output `d` contains the (50 x 4)-response matrix `d$Y`, 30 collinear predictors with an intercept stored in the matrix `d$X` as well as `d$R`, the simulated correlation between the responses.

**BVS4GCR** can be run using the command

    output <- BVS4GCR_MCMC(250, 50, 4, 100, d$Y, d$X, 2, c("Gaussian", "binary", "ordinal", "ordinal"), c(5, 3^2))

Inference with a 'dense' covariance matrix can be specified in the command line as follows

    output <- BVS4GCR_MCMC(250, 50, 4, 100, d$Y, d$X, 2, c("Gaussian", "binary", "ordinal", "ordinal"), c(5, 3^2), fullcov = TRUE)

## Example 3

This example coincides with the second simulated experiment presented in [Alexopoulos and Bottolo (2019)](https://arxiv.org/pdf/1907.08245.pdf) where 100 observations are simulated from a combination of one Gaussian, one binomial and two negative binomial distributions, with 100 predictors, an intercept and a common covariate for all responses

    set.seed(12345)
    d <- Sim_GCRdata(4, 100, 100, 2, c("Gaussian", "binomial", "negative binomial", "negative binomial"), extras = list(rhoY = 0.8, rhoX = 0.7, fixed_effect = c(-0.5, 0.5), mu_effect = c(1, rep(0.5, 3)), var_effect = c(1, rep(0.2, 3)), s1Lev = 0.05, s2Lev = 0.95, std_responses = rep(1, 4), binom_par = 10, negbin_par = c(0.5, 0.5)))

BVS is performed with **BVS4GCR** by typing

    output <- BVS4GCR_MCMC(250, 50, 4, 100, d$Y, d$X, 2, c("Gaussian", "binomial", "negative binomial", "negative binomial"), c(5, 3 ^2), ntrials = 10, negbin_init = c(0.5, 0.5))

The argument `nondecomp` allows to model the conditional dependence pattern between the responses with non-decomposable graphs 

    output <- BVS4GCR_MCMC(250, 50, 4, 100, d$Y, d$X, 2, c("Gaussian", "binomial", "negative binomial", "negative binomial"), c(5, 3 ^2), nondecomp = TRUE, ntrials = 10, negbin_init = c(0.5, 0.5))

## License and authors

This software uses the GPL v2 license, see [LICENSE](https://github.com/lb664/BVS4GCR/blob/master/LICENSE). Authors and copyright are provided in [DESCRIPTION](https://github.com/lb664/BVS4GCR/blob/master/DESCRIPTION)


## Issues

To report an issue, please use the **BVS4GCR** issue tracker at [BugReports](https://github.com/lb664/BVS4GCR/issues)

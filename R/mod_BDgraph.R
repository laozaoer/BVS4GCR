
# Internal function that performs one MCMC iteration to sample from the posterior distribution of a non-decomposable graph. The
# R code of this function is a modified version of the BDgraph R-package (Mohammadi and Wit, 2019). For details, see also Mohammadi 
# and Wit (2015) and Alexopoulos and Bottolo (2019)


########## mod_BDgraph ##########

mod_BDgraph <- function (data, n = NULL, niter = 5000, burnin = niter /2, not.cont = NULL, g.prior = 0.5, 
                         df.prior = 2, 
                         g.start, K.start, cores = NULL, threshold = 1e-08)
{
  if (niter < burnin)
  {
    stop("Number of iterations must be larger than the burn-in")
  }
  burnin <- floor(burnin)
  cores <- BDgraph::get_cores(cores = cores)
  list_S_n_p <- BDgraph::get_S_n_p(data = data, method = "ggm", n = n, not.cont = not.cont)
  S <- list_S_n_p$S
  n <- list_S_n_p$n
  p <- list_S_n_p$p
  method <- list_S_n_p$method
  colnames_data <- list_S_n_p$colnames_data
  b <- df.prior
  b_star <- b + n
  D <- diag(p)
  Ds <- D + S
  Ts <- chol(solve(Ds))
  Ti <- chol(solve(D))
  g_prior <- BDgraph::get_g_prior(g.prior = g.prior, p = p)
  G <- g.start
  G[g_prior == 1] <- 1
  G[g_prior == 0] <- 0
  G[lower.tri(G, diag(TRUE))] <- 0
  G <- G + t(G)
  K <- K.start
  p_links <- matrix(0, p, p)
  K_hat <- matrix(0, p, p)
  last_graph <- K_hat
  last_K <- K_hat
  
  result <- .C("ggm_bdmcmc_ma", as.integer(niter), as.integer(burnin),
               G = as.integer(G), as.double(g_prior), as.double(Ts),
               K = as.double(K), as.integer(p), as.double(threshold),
               K_hat = as.double(K_hat), p_links = as.double(p_links),
               as.integer(b), as.integer(b_star), as.double(Ds), as.integer(niter + 100), PACKAGE = "BDgraph")
  
  K_hat <- matrix(result$K_hat, p, p, dimnames = list(colnames_data, colnames_data))
  last_graph <- matrix(result$G, p, p, dimnames = list(colnames_data, colnames_data))
  last_K <- matrix(result$K, p, p)
  p_links = matrix(result$p_links, p, p, dimnames = list(colnames_data, colnames_data))
  p_links[lower.tri(p_links)] <- 0
  output <- list(p_links = p_links, last_K = last_K, last_graph = last_graph)
  
  return(output)
}


# Internal function for the proposal distribution of the binary latent vector for each response based on Guan and Stephens (2011). 
# For details, see Alexopoulos and Bottolo (2019)


########## Sample_gamma ##########

Sample_gamma <- function(cors, p, unoswap, gammacurr, uadd, geom_par)
{
  p_mix <- 0.3
  
  ########## comment ##########
  # Hard coding! p_mix is the mixture weight of the uniform component
  
  gam_idx <- which(gammacurr > 0)
  gam_zero_idx <- which(gammacurr == 0)
  pgam <- length(gam_idx)
  YXcorr_temp <- sort(cors, decreasing = TRUE, index.return = TRUE)
  YXcorr_rank_temp <- YXcorr_temp$ix
  YXcorr_rank <- intersect(YXcorr_rank_temp, gam_idx)
  YXcorr_zero_rank <- intersect(YXcorr_rank_temp, gam_zero_idx)
  gam_idx_prop <- gam_idx
  if ((runif(1) < unoswap) || (pgam == p) || (pgam == 0))
  {
    if (pgam > 0 & pgam < p)
    {
      if (runif(1) < uadd)
      {
        mymove <- "add"
        changej <- Qrnd(p, pgam, p_mix, geom_par)
        addj <- YXcorr_zero_rank[changej]
        gam_idx_prop <- union(gam_idx, addj)
        pgam_prop <- length(gam_idx_prop)
        lprop_gam <- (log(1 - uadd) - log(uadd)) + (-log(pgam_prop) - log(Qpdf(changej, p, pgam, p_mix, geom_par)))
      }else{
        mymove <- "delete"
        changej <- sample.int(pgam, 1)
        delj <- YXcorr_rank[changej]
        gam_idx_prop <- setdiff(gam_idx, delj)
        pgam_prop <- length(gam_idx_prop)
        lprop_gam <- (log(uadd) - log(1 - uadd)) + (log(pgam) + log(Qpdf(changej, p, pgam_prop, p_mix, geom_par)))
      }
    }
    if (pgam == 0 & pgam < p)
    {
      mymove <- "add"
      changej <- Qrnd(p, 0, p_mix, geom_par)
      addj <- YXcorr_zero_rank[changej]
      gam_idx_prop <- union(gam_idx, addj)
      pgam_prop <- length(gam_idx_prop)
      lprop_gam <- (log(1 - uadd) - log(1)) + (-log(pgam_prop) - log(Qpdf(changej, p, pgam, p_mix, geom_par)))
    }
    if (pgam > 0 & pgam == p){
      mymove <- "delete"
      changej <- sample.int(pgam, 1)
      delj <- YXcorr_rank[changej]
      gam_idx_prop <- setdiff(gam_idx, delj)
      pgam_prop <- length(gam_idx_prop)
      lprop_gam <- (log(uadd) -log(1) ) + (log(pgam) + log(Qpdf(changej, p, pgam_prop, p_mix, geom_par)))
      }
  }else{
    mymove <- "swap"
    changej <- Qrnd(p, pgam, p_mix, geom_par)
    addj <- YXcorr_zero_rank[changej]
    gam_idx_prop <- union(gam_idx, addj)
    pgam_prop <- length(gam_idx_prop)
    changej_OLD <- changej
    pgam_prop_OLD <- pgam_prop
    cond <- 0
    while (cond == 0)
    {
      changej <- sample.int(pgam_prop, 1)
      if (abs(changej - changej_OLD) > 0)
      {
        cond <- 1
      }
    }
    YX_corr_rank <- intersect(YXcorr_rank_temp, gam_idx_prop)
    delj <- YXcorr_rank[changej]
    gam_idx_prop <- setdiff(gam_idx, delj)
    pgam_prop <- length(gam_idx_prop)
    lprop_gam <- (-log(pgam_prop_OLD) - log(Qpdf(changej_OLD, p, pgam, p_mix, geom_par))) + (log(Qpdf(changej, p, pgam_prop, p_mix, geom_par)) + log(pgam_prop_OLD))
  }
  gam_idx_prop <- sort(gam_idx_prop)
  output <- list(gam_idx_prop, lprop_gam, mymove)
  
  return(output)
}

########## Qrnd ##########

Qrnd <- function(p, pgam, u, geom_par)
{
  if (runif(1) < u)
  {
    r <- sample.int(p - pgam, 1)
  }else{
    r <- qgeom(runif(1) * pgeom(p - pgam, 1 /geom_par), 1 /geom_par)
  }
  output <- r
  
  return(output)
}

########## Qpdf ##########

Qpdf <- function(r, p, pgam, u, geom_par)
{
  ff <- u * (1 /(p - pgam)) + (1 - u) * (dgeom(r, 1 /geom_par) / pgeom(p - pgam, 1 /geom_par))
  ouput <- ff
  
  return(ouput)
}

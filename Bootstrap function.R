##############################################
#Bootstrap function 
##############################################
bootstrap <- function(d)
{
  d[sample(x = seq(nrow(d)),
           size = nrow(d),
           replace = T),]
}

percentile_bootstrap <- function(d, A, estimator,ci_prob= c(0.025, 0.975))
{
  bs_data <- lapply(seq(A), function(x) bootstrap(d))
  bs_est <- lapply(bs_data, estimator)
  est_data <- sapply(bs_est, function(x) x[["psi10"]])
  ci <- quantile(est_data, probs = ci_prob)
  se <- sqrt(var(unlist( est_data))) 
  
  return(list("Bootstrap data" = bs_data,
              "lower_CI" = ci[1],
              "upper_CI"= ci[2],
              "SE" = se))
  
}


bootstrap_results  <- function (simulations, B, estimator, ci_prob= c(0.025, 0.975))
{
  bs_list  <- lapply(simulations,function(y) percentile_bootstrap(y, B, estimator,ci_prob))
  
  return(bs_list)
}

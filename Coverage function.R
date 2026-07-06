#########################################
#Coverage function
#########################################


get_coverage <- function(data_results, n_sims1,psi10_known)
{
  ci_low <- sapply(data_results, function(x) x["lower_CI"])
  ci_up <- sapply(data_results, function(x) x["upper_CI"])
  
  coverage<- NULL
  
  for(k in 1:n_sims1)
  {
    example_ci <- data.frame(ci_low[[k]],ci_up[[k]] )
    example_ci <- example_ci %>%
      mutate (cover = if_else ((example_ci[,1]< psi10_known &
                                  example_ci[,2]> psi10_known),1,0))
        coverage[k]<- example_ci$cover
      }
  
  coverage <- data.frame (coverage)
  
  cover_1 <- coverage %>% filter (coverage == 1)
  
  coverage_result <- nrow(cover_1)/nrow(coverage)
  
  return(coverage_result)
  
}
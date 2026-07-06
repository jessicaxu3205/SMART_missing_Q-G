######################################
#Q-learning function
######################################


hard_max_only <- function(data)
{
  data <- as.data.frame(data)
  
  data_r <- data[data$s==1,]
  
  q2   <- lm(Y ~ o1 + a1 + o1_a1 + a2 + o2_a2 + a1_a2, data = data_r)
  
  # Yprimeprime (also called the 'blip' function, or 'tailoring' vars)
  pseudo_Y_blip <- (coef(q2)[5]*data_r[["a2"]] +
                      # coef(q2)[6]*data_r[["o2_a2"]] +
                      coef(q2)[7]*data_r[["a1_a2"]])/data_r[["a2"]]

  # Yprime (also called the 'treatment free' function)
  pseudo_Y_free  <- coef(q2)[1] +
    coef(q2)[2]*data_r[["o1"]] +
    coef(q2)[3]*data_r[["a1"]] +
    coef(q2)[4]*data_r[["o1_a1"]]
  
  pseudo_Y    <- pseudo_Y_free +  abs(pseudo_Y_blip)
  data$y_s <- data$Y
  data$y_s[data$s==1] <- pseudo_Y
  q1          <- lm(y_s ~ o1 + a1 + o1_a1, data = data)
  
  return(q1)
}


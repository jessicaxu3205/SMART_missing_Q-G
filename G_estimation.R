#####################################################
#Fucntion for g_estimation
#####################################################

g_estimation <- function(data, conf= 0.95)
{
  
  data_t <- as.data.frame(data)
  
  n<- nrow(data_t)
  
  #Stage 2:------------------------------
  # treatment model 
  A2hat <- rep (0.5 , n )
  A2hat[data_t$s==0]<- 0
  
  # pseudo - outcome
  Y2 <- data_t$Y 
  # treatment - free model 
  H2.beta <- cbind ( rep (1 , n ) , data_t$o1, data_t$a1, data_t$o1_a1  )
  
  # blip model
  H2.psi <- cbind ( rep (1 , n ) , data_t$a1 )
  
  # w - matrix for analysis
  W2 <- diag ( data_t$a2 - A2hat )-( data_t$a2 - A2hat )*H2.beta %*% solve (t( H2.beta ) %*% H2.beta ) %*% t ( H2.beta )
  
  # estimate parameters
  psi2 <- solve ( t ( H2.psi ) %*% W2 %*% ( data_t$a2 * H2.psi ))%*% t ( H2.psi ) %*% W2 %*% Y2
  
  #Stage 1:------------------------------
  
  # treatment model 
  A1hat <- rep (0.5 , n )
  
  # pseudo - outcome
  Y1 <- Y2 - data_t$a2 *( H2.psi %*% psi2 )
  
  # treatment - free 
  H1.beta <- cbind ( rep (1 , n ),data_t$o1)
  
  # blip model
  H1.psi <- cbind ( rep (1 , n ), data_t$o1 )
  
  # weight - like matrix for analysis
  W1 <- diag ( data_t$a1 - A1hat ) - ( data_t$a1 - A1hat )* H1.beta%*% solve (t ( H1.beta ) %*% H1.beta ) %*% t ( H1.beta )
  
  # estimate parameters
  psi1 <- solve ( t ( H1.psi ) %*% W1 %*% ( data_t$a1 * H1.psi ))%*% t ( H1.psi ) %*% W1 %*% Y1
  
  #sandwich estimation------------------------------------------------------
  # stage 2:
  # bread = E[t(H2.psi)*A2*W2*H2.psi]
  B2 <- (1/n)*(t(H2.psi)%*%W2%*%(data_t$a2*H2.psi))
  
  # filling = E[U2*t(U2)]
  # where U2 is our estimating equation
  # for this, construct beta2 and diagonal matrix D2 with entries (A2-A2hat)
  beta2 <- solve(t(H2.beta) %*% H2.beta) %*% t(H2.beta) %*% (Y2 - data_t$a2*H2.psi%*%psi2)
  D2 <- diag(data_t$a2-A2hat)
  U2 <- t(D2 %*% H2.psi)%*%diag(as.vector(Y2 - data_t$a2*H2.psi %*% psi2 - H2.beta %*% beta2))
  F2 <- (1/n)*U2%*%t(U2)
  
  # sandwich covariance matrix = (B^-1) * F * t(B)^-1
  covmat2 <- (1/n) * solve(B2) %*% F2 %*% solve(t(B2))
  # standard error estimate for psi2 then the square root of the corresponding diagonal elements of
  #covmat2
  se2 <- sqrt(tail(diag(covmat2),2))
  
  # stage 1:
  # bread = E[t(H1.psi)*A1*W1*H1.psi]
  B1 <- (1/n)*(t(H1.psi)%*%W1%*%(data_t$a1*H1.psi))
  
  # filling = E[U1*t(U1)]
  # where U1 is our estimating equation
  # for this, construct beta1 and diagonal matrix D1 with entries (A1-A1hat)
  beta1 <- solve(t(H1.beta) %*% H1.beta) %*% t(H1.beta) %*% (Y1 - data_t$a1*H1.psi%*%psi1)
  D1 <- diag(data_t$a1-A1hat)
  U1 <- t(D1 %*% H1.psi)%*%diag(as.vector(Y1 - data_t$a1*H1.psi %*% psi1 - H1.beta %*% beta1))
  F1 <- (1/n)*U1%*%t(U1)
  
  # sandwich covariance matrix = (B^-1) * F * t(B)^-1
  covmat1 <- (1/n) * solve(B1) %*% F1 %*% solve(t(B1))
  # standard error estimate for psi1 then the square root of the corresponding diagonal elements of
  #covmat1
  
  se1 <- sqrt(tail(diag(covmat1),2))
  
  #confidence interval
  alpha<- (1-conf)/2
  z<-qnorm(1-alpha)
  CI<-psi1 +c(-1,1)*z*se1
  
  
  return(data.frame( "psi10"= psi1[1],
                     "sandwich_se"=se1[1],
                     "lower_CI" = CI[1],
                     "upper_CI"= CI[2]
  ))
  
}

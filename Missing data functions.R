#####################################################
#Missing data: functions used for multiple imputation
#####################################################
convert_fac<-function(data){
  
  data_m<- as.data.frame(data)
  
  data_m$a2 <- as.numeric(data_m$a2) 
  data_m$a2[data_m$a2==1] <- -1
  data_m$a2[data_m$a2==2] <- 0 
  data_m$a2[data_m$a2==3] <- 1 
  
  data_m$a1 <- as.numeric(data_m$a1) 
  data_m$a1[data_m$a1==1] <- -1
  data_m$a1[data_m$a1==2] <- 1  
  
  data_m$o2 <- as.numeric(data_m$o2) 
  data_m$o2[data_m$o2==1] <- -1
  data_m$o2[data_m$o2==2] <- 1  
  
  data_m$o1 <- as.numeric(data_m$o1) 
  data_m$o1[data_m$o1==1] <- -1
  data_m$o1[data_m$o1==2] <- 1  
  
  return(data_m)
}

#Pool estimates and se with Rubin's rule-----------------------------------
mi_pool <- function (mi_list, m, conf= 0.95)
{
  Q<- t(as.data.frame(lapply(mi_list,function(y) y[["psi10"]])))
  
  SE<-t(as.data.frame(lapply(mi_list,function(y) y[["sandwich_se"]])))
  
  #pooled estimate 
  Q_bar <-mean(Q)
  
  #Within-imputation variance
  U<- SE^2
  U_bar<- mean((U))
  
  #between variance
  B<- sum((Q-Q_bar)^2)/(m-1)
  
  #total variance 
  T<- U_bar +(1+1/m)*B
  
  SE_combined<- sqrt(T)
  
  #CI
  nu<- (m-1)*(1+(U_bar/((1+1/m)^B)))^2
  alpha<- (1-conf)/2
  t_crit<-qt(1-alpha,df=nu)
  CI<-Q_bar +c(-1,1)*t_crit*SE_combined
  
  return(data.frame ( pooled_est = Q_bar,
                      combined_SE = SE_combined,
                      lower_CI = CI[1],
                      upper_CI= CI[2]) 
  )
  
}

##Perform MI:----------------------------------------------

mi<- function(data, m,
              method)
{
  data_miss <- as.data.frame(data)
  
  data_m<- data_miss[,c(1:9)]
  
  data_m$o2_Y<-data_m$o2*data_m$Y
  data_m$a1_Y<-data_m$a1*data_m$Y
  data_m$a2_Y<-data_m$a2*data_m$Y
  
  
  data_m$o2 <- factor(data_m$o2) 
  data_m$o1 <- factor(data_m$o1) 
  data_m$a1 <- factor(data_m$a1) 
  
  data_m$a2 <- factor(data_m$a2) 
  
  pred <- make.predictorMatrix(data_m) 
  
  pred[,"o2_a2"]<- 0 ## o2_a2 is not used when imputing other variable, collinear
  pred[,"a2"]<- 0 ## a2 is not used when imputing other variable, collinear
  pred[,"o2"]<- 0 ## o2 is not used when imputing other variable, collinear
  
  pred["o2","o2_a2"]<- 0 ## o2_a2 is not used when imputing o2 variable
  pred["a2","o2_a2"]<- 0 ## o2_a2 is not used when imputing a2 variable
  
  pred["a2","a1_a2"]<- 0 ## a1_a2 is not used when imputing a2 variable
  
  pred[,"o2_Y"]<- 0 ## o2_Y is not used when imputing other variable
  pred[,"a1_Y"]<- 0 ## a1_Y is not used when imputing other variable
  pred[,"a2_Y"]<- 0 ## a2_Y is not used when imputing other variable
  
  pred["o2","a2_Y"]<- 1 ## a2_Y is used when imputing o2 variable
  pred["a2","o2_Y"]<- 1 ## o2_Y is used when imputing a2 variable
  pred["a2","a1_Y"]<- 1 ## a1_Y is used when imputing a2 variable
  
  
  mi_data <- data_m %>%
    mice( m=m, maxit=5,
          predictorMatrix=pred,
          method = method) %>% 
    mice::complete("all") %>%
    lapply(function(y) convert_fac(y))%>%
    lapply(function(y) g_estimation(y))
  
  pool <-mi_pool(mi_data,m)
  
  return(data.frame(pool))
  
}

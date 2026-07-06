##############################################
# Function to generate the simulated datasets
############################################


library(dplyr)
library(mice)
#####################################################
#Simulate dataset

#Set number of datasets
n_sim = 1000

#number of participants 
n=400

#Set up scenario 
gamma = c(0, 0, 0, 0, 0, 0, 0) #regular s1 -0.0004643346
delta = c(0.5,0.5)


expit <- function(x)
{
  exp(x)/(1 + exp(x))
}

###----------------------------------------------
Design1Simulate <- function(n, delta, n_sim){
  # Simulates n_sim design 1 SMART trials with n individuals with delta difference between best and non-best (confirm)
  #
  # Args: 
  #  n: Number of individuals in each SMART study simulation
  #  delta: Effect size
  #  n_sim: Number of simulations
  #
  # Returns: 
  # A list of n_sim SMART trial simulations each consisting of a data frame of the simulated variables, and outcomes
  ADHD.sim.list <- vector("list", length = n_sim)
  for (i in 1:n_sim){
    o1  <- rbinom(n, 1, 0.5)
    o1  <- 2*o1 - 1 # Converted to {-1,1} coding
    
    a1 <- 2*rbinom(n,1,0.5)-1 
    
    o2_prob <- expit(delta[1] * o1 +
                       delta[2] * a1)
    
    o2   <- rbinom(n, 1, o2_prob)
    
    s<-rep(0,n)
    s[o2==0]<-1 #Indicator of non-response
    
    a2<-2*rbinom(n,1,0.5)-1
    a2[s==0]<-0
    
    o2 <- 2*o2 - 1
    #####################
    
    o1_a1 <- o1*a1
    o2_a2 <- o2*a2
    a1_a2 <- a1*a2
    
    
    mu_y <- gamma[1] + gamma[2]*o1 + gamma[3]*a1 + gamma[4]*o1_a1 + 
      gamma[5]*a2 + gamma[6]*o2_a2 + gamma[7]*a1_a2
    
    # Observed vector.
    epsilon <- rnorm(n, 0, 1)
    Y <- mu_y + epsilon
    
    ADHD.sim.list[[i]]<-data.frame(o1,a1,o2,s,a2,o1_a1,o2_a2,a1_a2,Y)
    
  }
  ADHD.sim.list
}
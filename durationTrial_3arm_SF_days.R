rm(list=ls())

#### Allocation ratio for each treatment arm ####

library(dplyr)

# Durations of the treatment

durations <- c(6,4,3)

# Recruitment rate per month

recrate <- 30

# Total number of patients per arm needed for the trial

n <- 84*2

nc <- recrate/length(durations) #number of patients needed on the control 
ns <- nc #number of patients needed on the shortest arm 
nl <- nc #number of patients needed on the longest arm

# number of patients in each arm at the interim analysis

interim <- c(n/2+recrate/length(durations)*c(durations[1]-durations[2:3]),n/2)

# ratio between number of patients on each treatment arm and on the control

allocatioratio <- interim/interim[3]

# expected duration of the trial

expectdur <- ceiling(n/(recrate/length(durations)))

expectdur+durations[1]

#### calls functions (3-arm 2-stage design) with binary data

####  Function to compute boundary in ORD NI design for 3-arm and 2-stage #####

#' @u1: parameter for the grid search
#' @interimN: number of interim analyses
#' @alpha: alpha-level of the hypothesis test
#' @cov: covariance matrix
#' @first: covariance matrix at the first stage
#' @r: allocation ratio sample size in the first and second stage
#' @prec: precision for the boundaries: numbers after the comma
#' @arms: number of arms (3)
#' @ushape: upper boundary shape: "pocock", "obf", "triangular"
#' @lshape: lower boundary shape: "pocock", "obf", "triangular"

bounds_3arm2stageORD_binary <- function(u1, interimN, arms,r,
                                        alpha, cov, first, ushape,
                                        lshape,
                                        prec){
  
  library(mvtnorm)
  library(gtools)
  library(cubature)
  
  alphat <- rep(0,length = interimN)
  
  alphatcaseup<- NULL
  
  alphatcasedown<- NULL
  
  alphatcasefin<- NULL
  
  for(t in 1:length(u1)){
    
    for (i in 1:interimN){
      
      if (ushape == "obf") {
        u11 <- (u1[t] * 1/sqrt(r))[1:i]
      }
      else if (ushape == "pocock") {
        u11 <- rep(u1[t], i)
      }
      else if (ushape == "triangular") {
        u11 <- (u1[t] * (1 + r/max(r))/sqrt(r))[1:i]
      }
      
      if (lshape == "obf") {
        l11 <- c(-u1[t] * 1/sqrt(r))[1:(i-1)]
      }
      else if (lshape == "pocock") {
        l11 <- c(rep(-u1[t], i-1))
      }
      else if (lshape == "triangular") {
        if (ushape == "triangular") {
          l11 <- c(-u1[t] * ((1 -3* r/max(r))/sqrt(r)))[1:(i-1)]
        }
      }
      
      
      if(i ==1){
        set.seed(123456)
        
        alphat[i] <- pmvnorm(lower = u11[1],
                             upper = Inf,
                             sigma = first[1],
                             algorithm = GenzBretz(abseps = 1*10^-12))[1]
        
      }
      
      
      if(i>1){
        
        alphat[i] <- pmvnorm(lower = c(-Inf,u11[1], u11[2],-Inf),
                             upper = c(u11[1],Inf,Inf,Inf),
                             sigma = cov,
                             algorithm = GenzBretz(abseps = 1*10^-12))[1]+
          pmvnorm(lower = c(l11[1],l11[1], u11[2],-Inf),
                  upper = c(u11[1],u11[1],Inf,Inf),
                  sigma = cov,
                  algorithm = GenzBretz(abseps = 1*10^-12))[1]+
          pmvnorm(lower = c(l11[1],-Inf, u11[2],-Inf),
                  upper = c(u11[1],l11[1],Inf,Inf),
                  sigma = cov,
                  algorithm = GenzBretz(abseps = 1*10^-12))[1]
        
      }
      
    }
    
    
    if(sum(alphat) <= alpha){
      
      upperbound <- u11
      lowerbound <- l11
      
      totalpha <- sum(alphat)
      alow <- u1[t]
      finalalpha <- alphat
      
    }
    
    if(sum(alphat) > alpha){
      
      aup <- u1[t]
      break
      
      
    }
    
  }
  
  results <- list(upperbound, alow, aup, lowerbound,totalpha)
  
  names(results) <- c("upperbound", "alow", "aup", "lowerbound","totalpha")
  return(results)
  
}

#### Function to compute ss in ORD NI design for 3-arm and 2-stage ####

#' @theta: vector of clinically relevant difference 
#' @noninf: non-inferiority margin 
#' @stage: number of stages
#' @alpha: alpha-level of the hypothesis test
#' @beta: beta-level of the test
#' @r: allocation ratio sample size in the first and second stage
#' @rhonull: allocation ratio vector for the null configuration - (r_1^(1), r_1^(2), ..., r_1^(0))
#' @rhopower: allocation ratio vector for the power configuration - (r_1^(1), r_1^(2), ..., r_1^(0))
#' @prec: precision for the boundaries: numbers after the comma
#' @arms: number of arms (3)
#' @ushape: "pocock", "obf", "triangular"
#' @lshape: "pocock", "obf", "triangular"
#' @power: type of power - "reject at least one" or "reject all"

boundaries_3arms2stageORDNonInf_SF <- function(theta,
                                               noninf,
                                               stage,
                                               alpha, 
                                               beta,
                                               r, 
                                               rhonull1,
                                               prec,
                                               arms,
                                               ushape,
                                               lshape,
                                               power
){ 
  
  library(mvtnorm)
  library(cubature)
  
  pc <- theta[1]
  
  pt1 <- pc-noninf
  
  pt2 <- pc-noninf
  
  varpc <- pc*(1-pc)
  
  vart1 <- pt1*(1-pt1)
  
  vart2 <- pt2*(1-pt2)
  
  var <- c(vart1, vart2)
  
  pnull <- varpc+var
  
  rhonull <- rhopower <- c(1,1,1)
  
  covbfin12 <- (sqrt(rhonull1[1]*rhonull1[2])*pc*(1-pc))/(sqrt(rhonull1[3]*vart1+rhonull1[1]*varpc)*sqrt(rhonull1[3]*vart2+rhonull1[2]*varpc))
  
  covbfin13 <- sqrt(rhonull1[1]*rhonull1[3]*r[2]*rhonull[1]*r[2]*rhonull[3])/sqrt((r[2]*rhonull[3]*vart1+r[2]*rhonull[1]*varpc)*(r[1]*rhonull1[3]*vart1+r[1]*rhonull1[1]*varpc))*((vart1/(r[2]*rhonull[1]))+(varpc/(r[2]*rhonull[3])))
  
  covbfin14 <- sqrt(rhonull1[1]*rhonull1[3]*r[2]*rhonull[2]/(r[2]*rhonull[3]))*varpc/(sqrt(rhonull1[3]*vart1+rhonull1[1]*varpc)*sqrt(r[2]*rhonull[3]*vart2+r[2]*rhonull[2]*varpc))
  
  covbfin23 <- sqrt(rhonull1[2]*rhonull1[3]*r[2]*rhonull[1]/(r[2]*rhonull[3]))*varpc/(sqrt(rhonull1[3]*vart2+rhonull1[2]*varpc)*sqrt(r[2]*rhonull[3]*vart1+r[2]*rhonull[1]*varpc))
  
  covbfin24 <- sqrt(rhonull1[2]*rhonull1[3]*r[2]*rhonull[2]*r[2]*rhonull[3])/sqrt((r[2]*rhonull[3]*vart2+r[2]*rhonull[2]*varpc)*(r[1]*rhonull1[3]*vart2+r[1]*rhonull1[2]*varpc))*((vart2/(r[2]*rhonull[2]))+(varpc/(r[2]*rhonull[3])))
  
  covbfin34 <- (sqrt(r[2]*rhonull[1]*r[2]*rhonull[2])*pc*(1-pc))/(sqrt(r[2]*rhonull[3]*vart1+r[2]*rhonull[1]*varpc)*sqrt(r[2]*rhonull[3]*vart2+r[2]*rhonull[2]*varpc))
  
  covbfin <- matrix(c(1,covbfin12,covbfin13,covbfin14,
                      covbfin12, 1, covbfin23,covbfin24,
                      covbfin13,covbfin23,1,covbfin34,
                      covbfin14,covbfin24,covbfin34,1), nrow=4, byrow=TRUE)
  
  
  firststagematrix <- covbfin[1:(arms-1), 1:(arms-1)]
  
  u1 <- seq(from = 7, to = 0, by = -1)
  
  interimN <- stage
  
  first <- bounds_3arm2stageORD_binary(u1, interimN, arms,r, 
                                       alpha, covbfin, firststagematrix,
                                       ushape,lshape,prec)
  
  low <- first$alow
  
  up <- first$aup
  
  p <- rep(1, times = prec)
  
  for (p in 1:length(p)){
    
    callf <- bounds_3arm2stageORD_binary(seq(from = low, 
                                             to = up, by = -1/(10^p)), 
                                         interimN= interimN,
                                         arms,
                                         r,
                                         alpha,
                                         covbfin,
                                         firststagematrix,
                                         ushape,
                                         lshape,
                                         prec)
    
    low <- callf$alow
    
    up <- callf$aup
  }
  
  
  upperbound<- callf$upperbound
  
  lowerbound <- callf$lowerbound
  
  a <- callf$alow
  totalalpha <- callf$totalpha
  
  ########################## Power
  
  pc <- theta[1]
  
  pt1 <- theta[2]
  
  pt2 <- theta[3]
  
  varpc <- pc*(1-pc)
  
  vart1 <- pt1*(1-pt1)
  
  vart2 <- pt2*(1-pt2)
  
  var <- c(vart1, vart2)
  
  covbfin12 <- (sqrt(rhonull1[1]*rhonull1[2])*pc*(1-pc))/(sqrt(rhonull1[3]*vart1+rhonull1[1]*varpc)*sqrt(rhonull1[3]*vart2+rhonull1[2]*varpc))
  
  covbfin13 <- sqrt(rhonull1[1]*rhonull1[3]*r[2]*rhopower[1]*r[2]*rhopower[3])/sqrt((r[2]*rhopower[3]*vart1+r[2]*rhopower[1]*varpc)*(r[1]*rhonull1[3]*vart1+r[1]*rhonull1[1]*varpc))*((vart1/(r[2]*rhopower[1]))+(varpc/(r[2]*rhopower[3])))
  
  covbfin14 <- sqrt(rhonull1[1]*rhonull1[3]*r[2]*rhopower[2]/(r[2]*rhopower[3]))*varpc/(sqrt(rhonull1[3]*vart1+rhonull1[1]*varpc)*sqrt(r[2]*rhopower[3]*vart2+r[2]*rhopower[2]*varpc))
  
  covbfin23 <- sqrt(rhonull1[2]*rhonull1[3]*r[2]*rhopower[1]/(r[2]*rhopower[3]))*varpc/(sqrt(rhonull1[3]*vart2+rhonull1[2]*varpc)*sqrt(r[2]*rhopower[3]*vart1+r[2]*rhopower[1]*varpc))
  
  covbfin24 <- sqrt(rhonull1[2]*rhonull1[3]*r[2]*rhopower[2]*r[2]*rhopower[3])/sqrt((r[2]*rhopower[3]*vart2+r[2]*rhopower[2]*varpc)*(r[1]*rhonull1[3]*vart2+r[1]*rhonull1[2]*varpc))*((vart2/(r[2]*rhopower[2]))+(varpc/(r[2]*rhopower[3])))
  
  covbfin34 <- (sqrt(r[2]*rhopower[1]*r[2]*rhopower[2])*pc*(1-pc))/(sqrt(r[2]*rhopower[3]*vart1+r[2]*rhopower[1]*varpc)*sqrt(r[2]*rhopower[3]*vart2+r[2]*rhopower[2]*varpc))
  
  covbfinB <- matrix(c(1,covbfin12,covbfin13,covbfin14,
                       covbfin12, 1, covbfin23,covbfin24,
                       covbfin13,covbfin23,1,covbfin34,
                       covbfin14,covbfin24,covbfin34,1), nrow=4, byrow=TRUE)
  
  firststagematrixB <- covbfinB[1:(arms-1), 1:(arms-1)]
  
  maxsample <- NULL
  
  betat <- NULL
  
  betatcaseup <- NULL
  
  betatcasedown <- NULL
  
  betatcasefin <- NULL
  
  pop <- seq(from = 1, to = 4000, by = 1)
  
  for (n in 1:length(pop)){
    
    sigma1 <- (varpc/(rhonull1[3]*pop[n]))+(vart1/(rhonull1[1]*pop[n]))
    
    sigma2 <- (varpc/(rhonull1[3]*pop[n]))+(vart2/(rhonull1[2]*pop[n]))
    
    sigma21 <- (varpc/(rhopower[3]*r[2]*pop[n]))+(vart1/(rhopower[1]*r[2]*pop[n]))
    
    sigma22 <- (varpc/(rhopower[3]*r[2]*pop[n]))+(vart2/(rhopower[2]*r[2]*pop[n]))
    
    sigman <- c(sigma1,sigma2,sigma21,sigma22)
    
    mean2 <- rep(((theta[2:(arms)]-theta[1])+noninf),2)/(sqrt(sigman))
    
    mean1 <- mean2[1:(arms-1)]
    
    
    for (i in 1:interimN){
      
      u11 <- upperbound[1:i]
      l11 <- lowerbound[1:(i-1)]
      
      u2 <- upperbound[1:i]
      l2 <- lowerbound[1:(i-1)]
      
      if(power == "reject all"){
        
        if(i ==1){
          set.seed(123456)
          
          betat[i] <- pmvnorm(lower = c(u11[1],u2[1]),
                              upper = c(Inf, Inf),
                              sigma = firststagematrixB,
                              mean = mean1,
                              algorithm = GenzBretz(abseps = 1*10^-12))[1]
          
        }
        if(i>1){
          
          set.seed(123456)
          betat[i] <- pmvnorm(lower = c(l11[1],u2[1], u11[2],u2[2]),
                              upper = c(u11[1],Inf,Inf,Inf),
                              sigma = covbfinB,
                              mean = mean2,
                              algorithm = GenzBretz(abseps = 1*10^-12))[1]+
            pmvnorm(lower = c(-Inf,u2[1], u11[2],u2[2]),
                    upper = c(l11[1],Inf,Inf,Inf),
                    sigma = covbfinB,
                    mean = mean2,
                    algorithm = GenzBretz(abseps = 1*10^-12))[1]+
            pmvnorm(lower = c(u11[1],l2[1], -Inf,u2[2]),
                    upper = c(Inf,u2[1],Inf,Inf),
                    sigma = covbfinB,
                    mean = mean2,
                    algorithm = GenzBretz(abseps = 1*10^-12))[1]+
            pmvnorm(lower = c(l11[1],l2[1], u11[2],u2[2]),
                    upper = c(u11[1],u2[1],Inf,Inf),
                    sigma = covbfinB,
                    mean = mean2,
                    algorithm = GenzBretz(abseps = 1*10^-12))[1]
          
        }
      }
      if(power == "reject at least one"){
        
        if(i ==1){
          set.seed(123456)
          
          betat[i] <- pmvnorm(lower = u11[1],
                              upper = Inf,
                              sigma = firststagematrixB[1,1],
                              mean = mean1[1],
                              algorithm = GenzBretz(abseps = 1*10^-12))[1]
          
        }
        
        
        if(i>1){
          
          set.seed(123456)
          
          betat[i] <- pmvnorm(lower = c(-Inf,u11[1], u11[2],-Inf),
                              upper = c(u11[1],Inf,Inf,Inf),
                              sigma = covbfinB,
                              mean = mean2,
                              algorithm = GenzBretz(abseps = 1*10^-12))[1]+
            pmvnorm(lower = c(l11[1],l11[1], u11[2],-Inf),
                    upper = c(u11[1],u11[1],Inf,Inf),
                    sigma = covbfinB,
                    mean = mean2,
                    algorithm = GenzBretz(abseps = 1*10^-12))[1]+
            pmvnorm(lower = c(l11[1],-Inf, u11[2],-Inf),
                    upper = c(u11[1],l11[1],Inf,Inf),
                    sigma = covbfinB,
                    mean = mean2,
                    algorithm = GenzBretz(abseps = 1*10^-12))[1]
          
        }
        
        
        
      }
      
    }
    
    if(sum(betat) > 1-beta){
      
      maxsample <- pop[n]
      break
      
    }
  }
  
  results <- list(totalalpha,
                  lowerbound,upperbound,maxsample,a )
  
  names(results) <- c("sum_alphat", 
                      "lowerbounds", "upperbounds","sample size per arm per stage", "a")
  
  return(results)
}


#### SIMULATION 3-arm and 2-stage ORD NI design  for the SF strategy - match sample size at the final analysis ######

# PARAMETERS for both designs

seed <- 64736 # set seed to simulate the data

nsim <- 10^5 # number of simulations

arms <- 3 # number of arms

stage <- 2 # number of stages

alpha <- 0.05 # alpha-level of the test

noninf <- 0.1 # non-inferiority margin

pc <- 0.86 # response rate on the control arm

pt1 <- pc # response rate on the first treatment arm

pt2 <- pc # response rate on the second treatment arm

beta <- 0.2 # beta-level

ushape <- "triangular" # shape of upper critical boundary

lshape <- "triangular" # shape of lower critical boundary

thetapow <- 0 # value of treatment effect for power configuration

powerconf <- c(0, rep(thetapow, 2)) # power configuration

rhopower <- allocatioratio  #sample size in each treatment divided by sample size in the control: rho[1]=n11/n01, rho[2]=n21/n01, rho[3]=1

r <- ratio <- c(1,2) #allocation ratio sample size at different stages: n1=r[1], n2=r[2]*n1,...

prec <- 3 # precision - number of digits after decimal point - for the critical boundaries

# call function to find the sample size per arm per stage and the critical boundaries

fun <- boundaries_3arms2stageORDNonInf_SF(theta = c(pc,pt1,pt2)+powerconf,
                                          noninf,
                                          stage = stage,
                                          alpha = alpha,
                                          beta = beta,
                                          r = r,
                                          rhonull1=rhopower,
                                          prec = prec,
                                          arms = arms,
                                          ushape = ushape,
                                          lshape = ushape,
                                          power = "reject at least one")
fun

# probability to allocate every patient on each treatment arm

p <- c(10, 10, 10)/30  

### Scenario

### RESTRUCTURE scenarios

plac <- rep(pc, 2)

trt1 <- c(pc-noninf, pc)

trt2 <- c(pc-noninf, pc)

scen <- cbind(plac,trt1,trt2)

scen <- cbind(plac,trt1,trt2)

scenario <- matrix(data = scen,
                   nrow = nrow(scen),
                   ncol = ncol(scen))

scenario <- cbind(scenario, scenario[,2]-scenario[,1], scenario[,3]-scenario[,1])

colnames(scenario) <- c("p0", "p1", "p2","theta1", "theta2")

scenario <- as.data.frame(scenario)

######## Sample size and bound grid values 

pop <- fun$`sample size per arm per stage`#seq(fun$`sample size per arm per stage`-10, fun$`sample size per arm per stage`, by = 1)

upper <- fun$a#seq(fun$a-.2, fun$a+0.03, by = 0.01)

## simulations

upperfinal <- NULL
samplesize <- NULL

index_upper <- NULL
index_pop <- NULL
mat2 <- NULL

typeI <- power <- NULL

for(u in 1:length(pop)){
  
  popord <- pop[u]
  
  for(n in 1:length(upper)){
    
    #triangular bounds
    
    if(ushape == "triangular"){
      
      upp <- upper[n] * (1 + ratio/max(ratio))/sqrt(ratio)
      
      
      low <- -upper[n] * ((1 -3*ratio/max(ratio))/sqrt(ratio))[1]
    }
    
    #pocock bounds
    
    if(ushape == "pocock"){
      
      upp <- c(upper[n],upper[n])
      
      
      low <- -upper[n]
    }
    
    summary_ORD <- NULL
    
    u1 <- upp[1]
    
    u2 <- upp[2]
    
    l1 <- low[1]
    
    rej_all <- list()
    
    rej_all[[1]] <- matrix(data=FALSE, nrow=nrow(scenario),
                           ncol=nsim)
    
    rej_all[[2]] <- matrix(data=FALSE, nrow=nrow(scenario),
                           ncol=nsim)
    
    rej_any <- list()
    
    rej_any[[1]] <- matrix(data=FALSE, nrow=nrow(scenario),
                           ncol=nsim)
    
    rej_any[[2]] <- matrix(data=FALSE, nrow=nrow(scenario),
                           ncol=nsim)
    
    rej_ord_h01noh02 <- matrix(data=FALSE, nrow=nrow(scenario),
                               ncol=nsim)
    
    prop_ord <- list()
    
    prop_ord[[1]] <- matrix(data=NA, ncol=2,
                            nrow = nrow(scenario))
    prop_ord[[2]] <- matrix(data=NA, ncol=3,
                            nrow = nrow(scenario))
    
    estimatedss <- actualmaxss <-  matrix(data=NA, nrow=nrow(scenario),
                          ncol=nsim)
    
    estimatedsamplesize <- matrix(data=NA, ncol=1,
                                  nrow = nrow(scenario))
    
    durationtrial <- matrix(data=NA, ncol=1,
                            nrow = nrow(scenario))
    
    dur <- matrix(data=0, nrow=nrow(scenario),
                  ncol=nsim)
    
    durfin <- matrix(data=0, nrow=nrow(scenario),
                     ncol=nsim)
    
    
    for (k in 1:nrow(scenario)){
      
      set.seed(seed = seed)
      
      mean0 <- scenario[k,1]
      
      mean1 <- scenario[k,2]
      
      mean2 <- scenario[k,3]
      
      
      for (j in 1:nsim){
        
        plac1st <- dur1st <- dur2st <- NULL
        
        plac21 <- dur21 <- dur22 <- NULL
        
        plac_2missing <- plac <- 0
        
        count <- NULL
        
        #recruit until get popord on the placebo
        
        while(plac < popord){
          count <- c(count,sample(c("C","L","S"), 1, prob = p))
          plac <- sum(count=="C")
        }
        #table(count)
        
        n_plac <- sum(count=="C")
        n_long <- sum(count=="L")
        n_short <- sum(count=="S")
        days <- length(count)
        
        # continue to recruit until observe the patients on placebo
        
        count_over <- sample(c("C","L","S"), ceiling((durations[1]-durations[3])*30.41), 
                                          prob = p, replace=TRUE)
        
        #last_patient <- tapply(seq_along(count_over),count_over,max)
        
        n_plac_over <- sum(count_over=="C")
        n_long_over <- sum(count_over=="L")
        n_short_over <- sum(count_over=="S")
        
        n_long_overc <- sum(count_over[1:ceiling((durations[1]-durations[2])*30.41)]=="L")
        n_short_overc <- sum(count_over[1:ceiling((durations[1]-durations[3])*30.41)]=="S")
        
        plac11st <- rbinom(n = n_plac, prob= mean0, size = 1)
        
        dur11st <- rbinom(n = n_long+n_long_overc, prob= mean1, size = 1)
        
        dur21st <- rbinom(n = n_short+n_short_overc, prob= mean2, size = 1)
        
        dur[k,j] <- days + durations[1]
        
        estimatedss[k,j] <- length(plac11st)+length(dur11st)+length(dur21st)
        
        # Estimated response rates on the observed data
        
        p_plac11 <- mean(plac11st)
        
        p_dur11 <- mean(dur11st)
        
        p_dur12 <- mean(dur21st)
        
        
        p1 <- (mean1*(1-mean1)/length(dur11st))+(mean0*(1-mean0)/length(plac11st))
        
        p2 <- (mean2*(1-mean2)/length(dur21st))+(mean0*(1-mean0)/length(plac11st))
        
        
        sdpooled11 <- p1
        
        sdpooled12 <- p2
        
        # Z-statistics for the two treatment durations
        
        z1 <- (p_dur11-p_plac11+noninf)/sqrt(sdpooled11)
        
        z2 <- (p_dur12-p_plac11+noninf)/sqrt(sdpooled12)
        
        # The null hypothesis is early rejected in the 1st stage
        
        if( z1 >= u1 &  z2 >= u1){
          
          rej_all[[1]][k,j] = TRUE
          
        }
        
        if(z1 <= l1 &  z2 <= u1 ){
          
          rej_any[[1]][k,j] = TRUE
          
          
        }
        
        if(z1 >= u1 &  z2 <= l1){
          
          rej_ord_h01noh02[k,j] = TRUE
          
          
        }
        
        
        # Continue to second stage with all arms
        
        if((z1 <= u1 &   z2 >=u1)||
           ( z1 > l1  &  z1 < u1 &
             z2 > l1  &  z2 < u1)
        ){
          plac1st <- c(plac11st, rbinom(n = n_plac_over, prob= mean0, size = 1))
          
          dur1st <- c(dur11st, rbinom(n = n_long_over-n_long_overc, prob= mean1, size = 1))
          
          dur2st <- c(dur21st, rbinom(n = n_short_over-n_short_overc, prob= mean2, size = 1))
          
          # continue to recruit until observe the patients on the control duration
          
          plac_2missing <- 2*popord-(length(plac1st))
          
          plac_2 <- 0
          
          count_2 <- NULL
          
          #recruit until get popord on the placebo
          
          while(plac_2 < plac_2missing){
            count_2 <- c(count_2,sample(c("C","L","S"), 1, prob = p))
            days_2 <- length(count_2)
            plac_2 <- sum(count_2=="C")
          }
          #table(count_2)
          
          n_plac_2 <- sum(count_2=="C")
          n_long_2 <- sum(count_2=="L")
          n_short_2 <- sum(count_2=="S")
          
          # Simulate placebo and treatment populations
          
          plac21 <- c(plac1st, rbinom(n = n_plac_2, prob= mean0, size = 1))
          
          dur21 <- c(dur1st, rbinom(n = n_long_2, prob= mean1, size = 1))
          
          dur22 <- c(dur2st, rbinom(n = n_short_2, prob= mean2, size = 1))
          
          #duration
          
          dur[k,j] = dur[k,j]+days_2 + durations[1]
          
          #ESS
          
          estimatedss[k,j] <-  length(plac21)+length(dur21)+length(dur22)
          
          actualmaxss[k,j] <-  length(plac21)+length(dur21)+length(dur22)
          
          
          # Estimated response rates
          
          p_dur21 <- mean(dur21)
          
          p_plac21 <- mean(plac21)
          
          p_dur22 <- mean(dur22)
          
          
          p21 <- (mean1*(1-mean1)/length(dur21))+(mean0*(1-mean0)/length(plac21))
          
          p22 <- (mean2*(1-mean2)/length(dur22))+(mean0*(1-mean0)/length(plac21))
          
          
          sdpooled21 <- p21
          
          sdpooled22 <- p22
          
          
          # Z-statistics for the two treatment durations
          
          z21 <- (p_dur21-p_plac21+noninf)/sqrt(sdpooled21)
          
          z22 <- (p_dur22-p_plac21+noninf)/sqrt(sdpooled22)
          
          # Rejection of the null hypothesis at the second stage
          
          
          if( z21 >= u2 &  
              z22 >= u2){
            
            rej_all[[2]][k,j] = TRUE
            
          }
          
          if( z21 < u2){
            
            rej_any[[2]][k,j] = TRUE
            
            
          }
          
          
          if( z21 >= u2 &  
              z22 < u2){
            
            rej_ord_h01noh02[k,j] = TRUE
            
          }
          
        }
        
        # Continue to second stage with 1
        
        if((z1 > l1  &  z1 < u1 &
            z2 <  l1)
        ){
          
          plac1st <- c(plac11st, rbinom(n = n_plac_over, prob= mean0, size = 1))
          
          dur1st <- c(dur11st, rbinom(n = n_long_over-n_long_overc, prob= mean1, size = 1))
          
          # continue to recruit until observe the patients on the control duration
          
          plac_2missing <- 2*popord-(length(plac1st))
          
          plac_2 <- 0
          
          count_2 <- NULL
          
          #recruit until get popord on the placebo
          
          while(plac_2 < plac_2missing){
            count_2 <- c(count_2,sample(c("C","L"), 1, prob = p[1:2]))
            days_2 <- length(count_2)
            plac_2 <- sum(count_2=="C")
          }
          #table(count_2)
          
          n_plac_2 <- sum(count_2=="C")
          n_long_2 <- sum(count_2=="L")
          
          # Simulate placebo and treatment populations
          
          plac21 <- c(plac1st, rbinom(n = n_plac_2, prob= mean0, size = 1))
          
          dur21 <- c(dur1st, rbinom(n = n_long_2, prob= mean1, size = 1))
          
          #duration
          
          dur[k,j] = dur[k,j]+days_2 + durations[1]
          
          #ESS
          
          estimatedss[k,j] <- length(plac21)+length(dur21)
          
          
          # Estimated response rates
          
          p_plac21 <- mean(plac21)
          
          p_dur21 <- mean(dur21)
          
          p21 <- (mean1*(1-mean1)/length(dur21))+(mean0*(1-mean0)/length(plac21))
          
          sdpooled21 <- p21
          
          # Z-statistics for the two treatment durations
          
          z21 <- (p_dur21-p_plac21+noninf)/sqrt(sdpooled21)
          
          # Rejection of the null hypothesis at the second stage
          
          if (z21 < u2){
            
            rej_any[[2]][k,j] = TRUE
            
          }
          
          if( z21 >= u2){
            
            rej_ord_h01noh02[k,j] = TRUE
            
          }
          
        }
        
        # Continue to second stage with 2
        
        if((z1 >= u1 &
            z2 <  u1 &  z2 > l1)
        ){
          
          plac1st <- c(plac11st, rbinom(n = n_plac_over, prob= mean0, size = 1))
          
          dur2st <- c(dur21st, rbinom(n = n_short_over-n_short_overc, prob= mean2, size = 1))
          
          # continue to recruit until observe the patients on the control duration
          
          plac_2missing <- 2*popord-(length(plac1st))
          
          plac_2 <- 0
          
          count_2 <- NULL
          
          #recruit until get popord on the placebo
          
          while(plac_2 < plac_2missing){
            count_2 <- c(count_2,sample(c("C","S"), 1, prob = p[c(1,3)]))
            days_2 <- length(count_2)
            plac_2 <- sum(count_2=="C")
          }
          #table(count_2)
          
          n_plac_2 <- sum(count_2=="C")
          
          n_short_2 <- sum(count_2=="S")
          
          # Simulate placebo and treatment populations
          
          plac21 <- c(plac1st, rbinom(n = n_plac_2, prob= mean0, size = 1))
          
          dur22 <- c(dur2st, rbinom(n = n_short_2, prob= mean2, size = 1))
          
          #duration
          
          dur[k,j] = dur[k,j]+days_2 + durations[1]
          
          #ESS
          estimatedss[k,j] <- length(plac21)+length(dur22)
          
          # Estimated response rates
          
          p_plac21 <- mean(plac21)
          
          p_dur22 <- mean(dur22)
          
          p22 <- (mean2*(1-mean2)/length(dur22))+(mean0*(1-mean0)/length(plac21))
          
          
          sdpooled22 <- p22
          
          # Z-statistics for the two treatment durations
          
          z22 <- (p_dur22-p_plac21+noninf)/sqrt(sdpooled22)
          
          # Rejection of the null hypothesis at the second stage
          
          
          if( z22 >= u2){
            
            rej_all[[2]][k,j] = TRUE
            
          }
          
          if( z22 < u2){
            
            rej_ord_h01noh02[k,j] = TRUE
            
          }
          
        }
        
        
      }
      # Average probabilities of rej the hypotheses
      
      prop_ord[[1]][k,1] <- sum(rej_all[[1]][k,])/nsim
      
      prop_ord[[1]][k,2] <- sum(rej_any[[1]][k,])/nsim
      
      prop_ord[[2]][k,1] <- sum(rej_all[[2]][k,])/nsim
      
      prop_ord[[2]][k,2] <- sum(rej_any[[2]][k,])/nsim
      
      prop_ord[[2]][k,3] <- sum(rej_ord_h01noh02[k,])/nsim
      
      estimatedsamplesize[k] <- mean(estimatedss[k,])
      
      durationtrial[k] <- mean(dur[k,])
      
    }
    
    
    
    colnames(prop_ord[[1]]) <- c("ORD_rejall_1stage",
                                 "ORD_rejany_1stage")
    
    colnames(prop_ord[[2]]) <- c("ORD_rejall_2stage",
                                 "ORD_rejany_2stage",
                                 "ORD_rejH01noH02")
    
    colnames(estimatedsamplesize) <- "ESS Ordered ORD"
    
    colnames(durationtrial) <- "Expected duration (days)"
    
    summary_ORD <- cbind(popord,
                         2*popord*arms,
                         c(mean(actualmaxss[1,], na.rm= TRUE),mean(actualmaxss[2,], na.rm= TRUE)),#ceiling(popord*2*arms),
                         scenario,
                         prop_ord,
                         estimatedsamplesize,
                         durationtrial,
                         u1,
                         u2,
                         l1, 
                         allocatioratio[1],allocatioratio[2])
    
    colnames(summary_ORD)[1:3]<- c("patientscontrol_perstage ORD",
                                   "Max theoretical SS",
                                   "Actual SS")
    
    colnames(summary_ORD)[19:20]<- c("alloc ratio T_L/T_c",
                                     "alloc ratio T_S/T_c")
    
    summary_ORD <- summary_ORD %>% mutate(
      
      rej_all_overall = ORD_rejall_1stage + ORD_rejall_2stage,
      
      rej_any_overall = ORD_rejany_1stage + ORD_rejany_2stage,
      
      typeI = ORD_rejH01noH02 + rej_all_overall,
      
      sumprob = typeI + rej_any_overall,
      
      `Expected duration (months)` = durationtrial/30.41
      
    )%>% mutate(
      
      prob = paste("(", p0, ", ", p1,", ",p2, ")", sep = ""),
      theta = paste("(", round(theta1, 3), ", ", round(theta2, 3), ")", sep = "")
      
    )%>% select(
      
      `patientscontrol_perstage ORD`,
      `Max theoretical SS`,
      `Actual SS`,
      `ESS Ordered ORD`,
      prob,
      theta,
      typeI,
      rej_all_overall,
      rej_any_overall,
      `Expected duration (days)`,
      `Expected duration (months)`,
      sumprob,
      u1,
      u2,
      l1,
      `alloc ratio T_L/T_c`,
      `alloc ratio T_S/T_c`
      
      
    )
    
    colnames(summary_ORD) <- c("Patients in control arm (1st stage)", 
                               "Max theoretical SS", 
                               "Actual Max SS",
                               "ESS",
                               "Response rate (p0, p1, p2)",
                               "Treatment effect (theta1, theta2)",
                               "Reject at least one hyp.",
                               "Reject all hyp.",
                               "Reject neither hyp.",
                               "Expected duration (days)",
                               "Expected duration (months)",
                               "sumprob",
                               "u1",
                               "u2",
                               "l1",
                               "alloc ratio T_L/T_C",
                               "alloc ratio T_S/T_C")
    
    typeI <- summary_ORD[1,]$`Reject at least one hyp.`
    
    power <- summary_ORD[2,]$`Reject at least one hyp.`
    
    if(typeI <= alpha & power > 1-beta ){
      
      mat2 <- rbind(mat2, cbind(pop[u],
                                upper[n],
                                u1,
                                u2,
                                l1,
                                typeI,
                                power))
      
      
      
    }
    
    
    
  }
  
  
}

colnames(mat2) <- c("sample size", "a","u1", "u2", "l1", "typeI under null", "power")

mat2 <- cbind(mat2, pc,noninf)


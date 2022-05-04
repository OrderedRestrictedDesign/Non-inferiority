rm(list=ls())

library(dplyr)

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
                             sigma = first[1])[1]
        
      }
      
      
      if(i>1){
        
        alphat[i] <- pmvnorm(lower = c(-Inf,u11[1], u11[2],-Inf),
                             upper = c(u11[1],Inf,Inf,Inf),
                             sigma = cov)[1]+
          pmvnorm(lower = c(l11[1],l11[1], u11[2],-Inf),
                  upper = c(u11[1],u11[1],Inf,Inf),
                  sigma = cov)[1]+
          pmvnorm(lower = c(l11[1],-Inf, u11[2],-Inf),
                  upper = c(u11[1],l11[1],Inf,Inf),
                  sigma = cov)[1]
        
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

boundaries_3arms2stageORDNonInf <- function(theta,
                                            noninf,
                                            stage,
                                            alpha, 
                                            beta,
                                            r, 
                                            rhonull,
                                            rhopower, 
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
  
  covbfin12 <- (sqrt(rhonull[1]*rhonull[2])*pc*(1-pc))/(sqrt(rhonull[3]*vart1+rhonull[1]*varpc)*sqrt(rhonull[3]*vart2+rhonull[2]*varpc))
  
  covbfin13 <- sqrt(rhonull[1]*rhonull[3]*r[2]*rhonull[1]*r[2]*rhonull[3])/sqrt((r[2]*rhonull[3]*vart1+r[2]*rhonull[1]*varpc)*(r[1]*rhonull[3]*vart1+r[1]*rhonull[1]*varpc))*((vart1/(r[2]*rhonull[1]))+(varpc/(r[2]*rhonull[3])))
  
  covbfin14 <- sqrt(rhonull[1]*rhonull[3]*r[2]*rhonull[2]/(r[2]*rhonull[3]))*varpc/(sqrt(rhonull[3]*vart1+rhonull[1]*varpc)*sqrt(r[2]*rhonull[3]*vart2+r[2]*rhonull[2]*varpc))
  
  covbfin23 <- sqrt(rhonull[2]*rhonull[3]*r[2]*rhonull[1]/(r[2]*rhonull[3]))*varpc/(sqrt(rhonull[3]*vart2+rhonull[2]*varpc)*sqrt(r[2]*rhonull[3]*vart1+r[2]*rhonull[1]*varpc))
  
  covbfin24 <- sqrt(rhonull[2]*rhonull[3]*r[2]*rhonull[2]*r[2]*rhonull[3])/sqrt((r[2]*rhonull[3]*vart2+r[2]*rhonull[2]*varpc)*(r[1]*rhonull[3]*vart2+r[1]*rhonull[2]*varpc))*((vart2/(r[2]*rhonull[2]))+(varpc/(r[2]*rhonull[3])))
  
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
  
  covbfin12 <- (sqrt(rhopower[1]*rhopower[2])*pc*(1-pc))/(sqrt(rhopower[3]*vart1+rhopower[1]*varpc)*sqrt(rhopower[3]*vart2+rhopower[2]*varpc))
  
  covbfin13 <- sqrt(rhopower[1]*rhopower[3]*r[2]*rhopower[1]*r[2]*rhopower[3])/sqrt((r[2]*rhopower[3]*vart1+r[2]*rhopower[1]*varpc)*(r[1]*rhopower[3]*vart1+r[1]*rhopower[1]*varpc))*((vart1/(r[2]*rhopower[1]))+(varpc/(r[2]*rhopower[3])))
  
  covbfin14 <- sqrt(rhopower[1]*rhopower[3]*r[2]*rhopower[2]/(r[2]*rhopower[3]))*varpc/(sqrt(rhopower[3]*vart1+rhopower[1]*varpc)*sqrt(r[2]*rhopower[3]*vart2+r[2]*rhopower[2]*varpc))
  
  covbfin23 <- sqrt(rhopower[2]*rhopower[3]*r[2]*rhopower[1]/(r[2]*rhopower[3]))*varpc/(sqrt(rhopower[3]*vart2+rhopower[2]*varpc)*sqrt(r[2]*rhopower[3]*vart1+r[2]*rhopower[1]*varpc))
  
  covbfin24 <- sqrt(rhopower[2]*rhopower[3]*r[2]*rhopower[2]*r[2]*rhopower[3])/sqrt((r[2]*rhopower[3]*vart2+r[2]*rhopower[2]*varpc)*(r[1]*rhopower[3]*vart2+r[1]*rhopower[2]*varpc))*((vart2/(r[2]*rhopower[2]))+(varpc/(r[2]*rhopower[3])))
  
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
    
    sigma1 <- (varpc/(rhopower[3]*pop[n]))+(vart1/(rhopower[1]*pop[n]))
    
    sigma2 <- (varpc/(rhopower[3]*pop[n]))+(vart2/(rhopower[2]*pop[n]))
    
    sigman <- c(sigma1,sigma2)*rep(1/r, each=arms-1)
    
    mean2<- rep(((theta[2:(arms)]-theta[1])+noninf),2)/(sqrt(sigman))
    
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
                              mean = mean1)[1]
          
        }
        if(i>1){
          
          set.seed(123456)
          betat[i] <- pmvnorm(lower = c(l11[1],u2[1], u11[2],u2[2]),
                              upper = c(u11[1],Inf,Inf,Inf),
                              sigma = covbfinB,
                              mean = mean2)[1]+
            pmvnorm(lower = c(-Inf,u2[1], u11[2],u2[2]),
                    upper = c(l11[1],Inf,Inf,Inf),
                    sigma = covbfinB,
                    mean = mean2)[1]+
            pmvnorm(lower = c(u11[1],l2[1], -Inf,u2[2]),
                    upper = c(Inf,u2[1],Inf,Inf),
                    sigma = covbfinB,
                    mean = mean2)[1]+
            pmvnorm(lower = c(l11[1],l2[1], u11[2],u2[2]),
                    upper = c(u11[1],u2[1],Inf,Inf),
                    sigma = covbfinB,
                    mean = mean2)[1]
          
        }
      }
      if(power == "reject at least one"){
        
        if(i ==1){
          set.seed(123456)
          
          betat[i] <- pmvnorm(lower = u11[1],
                              upper = Inf,
                              sigma = firststagematrixB[1,1],
                              mean = mean1[1])[1]
          
        }
        
        
        if(i>1){
          
          set.seed(123456)
          betat[i] <- pmvnorm(lower = c(l11[1],-Inf, u11[2],-Inf),
                              upper = c(u11[1],Inf,Inf,Inf),
                              sigma = covbfinB,
                              mean = mean2)[1]+
            pmvnorm(lower = c(-Inf,u2[1], u11[2],-Inf),
                    upper = c(l11[1],Inf,Inf,Inf),
                    sigma = covbfinB,
                    mean = mean2)[1]
          
          
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

#### SIMULATION 3-arm and 2-stage ORD NI design ######

# PARAMETERS for both designs

seed <- 64736 # set seed to simulate the data

nsim <- 10^6 # number of simulations

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

rho <- c(1,1,1) #sample size in each treatment divided by sample size in the control: rho[1]=n11/n01, rho[2]=n21/n01, rho[3]=1

r <- ratio <- c(1,2) #allocation ratio sample size at different stages: n1=r[1], n2=r[2]*n1,...

prec <- 3 # precision - number of digits after decimal point - for the critical boundaries

# call function to find the sample size per arm per stage and the critical boundaries

fun <- boundaries_3arms2stageORDNonInf(theta = c(pc,pt1,pt2)+powerconf,
                                       noninf,
                                       stage = stage,
                                       alpha = alpha,
                                       beta = beta,
                                       r = r,
                                       rhonull=rho,
                                       rhopower=rho, 
                                       prec = prec,
                                       arms = arms,
                                       ushape = ushape,
                                       lshape = ushape,
                                       power = "reject at least one")
fun

### Scenario

### RESTRUCTURE scenarios

plac <- rep(pc, 2)

trt1 <- c(pc-noninf, pc)

trt2 <-  c(pc-noninf, pc)

scen <- cbind(plac,trt1,trt2)

scenario <- matrix(data = scen,
                   nrow = nrow(scen),
                   ncol = ncol(scen))

scenario <- cbind(scenario, scenario[,2]-scenario[,1], scenario[,3]-scenario[,1])

colnames(scenario) <- c("p0", "pL", "pS","thetaL", "thetaS")

scenario <- as.data.frame(scenario)

######## Sample size and bound grid values 

upper <- fun$a #seq(round(fun$alow,3)-0.1,round(fun$alow,3)+0.1,by = 0.001)

pop <- fun$`sample size per arm per stage`#seq(fun$`sample size per arm per stage`-5,fun$`sample size per arm per stage`+3,by=1)

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
    
    # critical boundaries
    
    u1 <- ifelse(ushape == "obf", (upper[n] * 1/sqrt(ratio))[1],
                 ifelse(ushape == "pocock",rep(upper[n], stage-1),
                        (upper[n] * (1 + ratio/max(ratio))/sqrt(ratio))[1]))
    
    l1 <- ifelse(ushape == "obf", (-upper[n] * 1/sqrt(ratio)),
                 ifelse(ushape == "pocock",rep(-upper[n], stage-1),
                        (-upper[n] * ((1 -3* ratio/max(ratio))/sqrt(ratio)))[1]))
    
    u2 <- ifelse(ushape == "obf", (upper[n] * 1/sqrt(ratio))[2],
                 ifelse(ushape == "pocock",rep(upper[n], stage-1),
                        (upper[n] * (1 + ratio/max(ratio))/sqrt(ratio))[2]))
    
    summary_ORD <- NULL
    
    # Save when rejecting any hypotheses

    bothineff <- matrix(data=FALSE, nrow=nrow(scenario),
                        ncol=nsim)
    
    # Save rejections of rejecting H01 and H02
    
    rej_ord_h01andh02 <- list()
    
    rej_ord_h01andh02[[1]] <- matrix(data=FALSE, nrow=nrow(scenario),
                                     ncol=nsim)
    
    rej_ord_h01andh02[[2]] <- matrix(data=FALSE, nrow=nrow(scenario),
                                     ncol=nsim)
    
    # Save rejections of rejecting H01 and not H02
    
    rej_ord_h01noh02 <- list()
    
    rej_ord_h01noh02[[1]] <- matrix(data=FALSE, nrow=nrow(scenario),
                                    ncol=nsim)
    
    rej_ord_h01noh02[[2]] <- matrix(data=FALSE, nrow=nrow(scenario),
                                    ncol=nsim)
    
    # Save all rejections for each scenario
    
    prop_ord <- list()
    
    prop_ord[[1]] <- matrix(data=NA, ncol=2,
                            nrow = nrow(scenario))
    prop_ord[[2]] <- matrix(data=NA, ncol=3,
                            nrow = nrow(scenario))
    
    # Estimated sample size for each scenario. --> The average of the sample size used for each simulation
    
    estimatedss <- matrix(data=NA, nrow=nrow(scenario),
                          ncol=nsim)
    
    estimatedsamplesize <- matrix(data=NA, ncol=1,
                                  nrow = nrow(scenario))

for (k in 1:nrow(scenario)){

  set.seed(seed = seed)

  mean0 <- scenario[k,1]

  mean1 <- scenario[k,2]

  mean2 <- scenario[k,3]

  for (j in 1:nsim){

    plac1st <- rbinom(n = popord, size = 1, prob = mean0)

    dur1st <- rbinom(n = popord, size = 1, prob = mean1)

    dur2st <- rbinom(n = popord, size = 1, prob = mean2)

    estimatedss[k,j] <- length(plac1st)+length(dur1st)+length(dur2st)


    # Estimated response rates on the observed data

    p_plac11 <- mean(plac1st)

    p_dur11 <- mean(dur1st)

    p_dur12 <- mean(dur2st)
    
    var1 <- (mean1*(1-mean1)/length(dur1st))+(mean0*(1-mean0)/length(plac1st))
    
    var2 <- (mean2*(1-mean2)/length(dur2st))+(mean0*(1-mean0)/length(plac1st))

    # Z-statistics for the two treatment durations

    z1 <- (p_dur11-p_plac11+noninf)/sqrt(var1)

    z2 <- (p_dur12-p_plac11+noninf)/sqrt(var2)

    # The null hypothesis is early rejected in the 1st stage

  
    if( z1 >= u1 &  z2 >= u1){

      rej_ord_h01andh02[[1]][k,j] = TRUE

    }


    if( z1 >= u1 &  z2 < l1){

      rej_ord_h01noh02[[1]][k,j] = TRUE


    }

    if( z1 < l1 &  z2 < u1){

      bothineff[k,j] = TRUE


    }

    # Continue to second stage with both treatment arms
    
    if((z1 < u1 & z1 >= l1 &
       z2 < u1 & z2 >= l1) ||
       (z1 < u1 &
        z2 >= u1)){

      # Simulate placebo and treatment populations

      plac21_2 <- rbinom(n = popord, size = 1, prob = mean0)

      dur21_2 <- rbinom(n = popord, size = 1, prob = mean1)

      dur22_2 <- rbinom(n = popord, size = 1, prob = mean2)

      estimatedss[k,j] <-  estimatedss[k,j]+length(plac21_2)+length(dur21_2)+length(dur22_2)

      plac21 <- c(plac21_2, plac1st)

      dur21 <- c(dur21_2, dur1st)

      dur22 <- c(dur22_2, dur2st)

      # Estimated response rates

      p_dur21 <- mean(dur21)

      p_plac21 <- mean(plac21)

      p_dur22 <- mean(dur22)
      
      var21 <- (mean1*(1-mean1)/length(dur21))+(mean0*(1-mean0)/length(plac21))
      
      var22 <- (mean2*(1-mean2)/length(dur22))+(mean0*(1-mean0)/length(plac21))
      
      # Z-statistics for the two treatment durations
      
      z21 <- (p_dur21-p_plac21+noninf)/sqrt(var21)

      z22 <- (p_dur22-p_plac21+noninf)/sqrt(var22)

      # Rejection of the null hypothesis at the second stage

      if (z21 >=  u2 &  z22 >= u2){

        rej_ord_h01andh02[[2]][k,j] = TRUE


      }

      if (z21 >= u2  & z22 < u2){

        rej_ord_h01noh02[[2]][k,j] = TRUE


      }


      if( z21 < u2){

        bothineff[k,j] = TRUE


      }

    }

    # Continue only with second treatment arm

    if(z1 >= u1 &
       z2 < u1 & z2 > l1){

      # Simulate placebo and treatment populations

      plac21_2 <- rbinom(n = popord, size = 1, prob = mean0)

      dur22_2 <- rbinom(n = popord, size = 1, prob = mean2)

      estimatedss[k,j] <- estimatedss[k,j]+length(plac21_2)+length(dur22_2)

      plac21 <- c(plac21_2, plac1st)

      dur22 <- c(dur22_2, dur2st)

      # Estimated response rates

      p_plac21 <- mean(plac21)

      p_dur22 <- mean(dur22)

      var22 <- (mean2*(1-mean2)/length(dur22))+(mean0*(1-mean0)/length(plac21))
      
      # Z-statistics for the two treatment durations
      
      z22 <- (p_dur22-p_plac21+noninf)/sqrt(var22)

      # Rejection of the null hypothesis at the second stage

      if (z22 >=  u2){

       
        rej_ord_h01andh02[[2]][k,j] = TRUE


      }

      if (z22 < u2){

        rej_ord_h01noh02[[2]][k,j] = TRUE

      }


    }

    # Continue only with first treatment arm

    if(z2 <= l1 &
       z1 < u1 & z1 > l1){

      plac21_2 <- rbinom(n = popord, size = 1, prob = mean0)

      dur21_2 <- rbinom(n = popord, size = 1, prob = mean1)

      estimatedss[k,j] <- estimatedss[k,j]+length(plac21_2)+length(dur21_2)

      plac21 <- c(plac21_2, plac1st)

      dur21 <- c(dur21_2, dur1st)

      # Estimated response rates

      p_plac21 <- mean(plac21)

      p_dur21 <- mean(dur21)

      var21 <- (mean1*(1-mean1)/length(dur21))+(mean0*(1-mean0)/length(plac21))
      
      # Z-statistics for the two treatment durations
      
      z21 <- (p_dur21-p_plac21+noninf)/sqrt(var21)
      

      # Rejection of the null hypothesis at the second stage

      if (z21 >= u2){

        rej_ord_h01noh02[[2]][k,j] = TRUE

      }
      if( z21 < u2){

        bothineff[k,j] = TRUE


      }


    }


  }
  # Average probabilities of rej the hypotheses

  prop_ord[[1]][k,1] <- sum(rej_ord_h01andh02[[1]][k,])/nsim

  prop_ord[[1]][k,2] <- sum(rej_ord_h01noh02[[1]][k,])/nsim

  prop_ord[[2]][k,1] <- sum(rej_ord_h01andh02[[2]][k,])/nsim

  prop_ord[[2]][k,2] <- sum(rej_ord_h01noh02[[2]][k,])/nsim

  prop_ord[[2]][k,3] <- sum(bothineff[k,])/nsim

  estimatedsamplesize[k] <- mean(estimatedss[k,])

}

    colnames(prop_ord[[1]]) <- c(
                                 "ORD_rejH01andH02_1stage",
                                 "ORD_rejH01noH02_1stage")
    
    colnames(prop_ord[[2]]) <- c(
                                 "ORD_rejH01andH02_2stage",
                                 "ORD_rejH01noH02_2stage",
                                
                                 "ORD_rejneitherH01andH02"
                                )
    
    colnames(estimatedsamplesize) <- "ESS Ordered ORD"
    
    summary_ORD <- cbind(popord,
                         ceiling(popord*2*arms),
                         scenario,
                         prop_ord,
                         estimatedsamplesize,
                         u1,
                         u2,
                         l1) %>% mutate(
                           
                           Ord_rejall = ORD_rejH01andH02_1stage+ORD_rejH01andH02_2stage,
                           Ord_rej1not2 = ORD_rejH01noH02_1stage+ORD_rejH01noH02_2stage,
                           
                           typeI= Ord_rejall+Ord_rej1not2,
                           sumprob = typeI+ORD_rejneitherH01andH02
                         )
    
    colnames(summary_ORD)[1:2]<- c("patients_perarm1stage ORD",
                                   "total sample ORD")
    
    typeI <- summary_ORD[1,]$typeI
    
    power <- summary_ORD[2,]$typeI
    
    if(typeI <= alpha & power > 1-beta){
      
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

mat2

summary_ORD
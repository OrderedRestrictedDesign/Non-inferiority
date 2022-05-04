rm(list=ls())

library(dplyr)

### Function to find critical boundaries for 4-arm 2-stage ORD NI design ####

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

bounds_4arm2stageORD <- function(u, arms,ratio,
                                 alpha, cov, first, ushape,lshape,prec){
  
  library(mvtnorm)
  library(gtools)
  library(cubature)
  
  for(z in 1:length(u)){
    
    prob <- 0
    
    if (ushape == "obf") {
      u1 <- (u[z] * 1/sqrt(ratio))[1]
    }
    else if (ushape == "pocock") {
      u1 <- rep(u[z], 1)
    }
    else if (ushape == "triangular") {
      u1 <- (u[z] * (1 + ratio/max(ratio))/sqrt(ratio))[1]
    }
    
    if (lshape == "obf") {
      l1 <- c(-u[z] * 1/sqrt(ratio))[1]
    }
    else if (lshape == "pocock") {
      l1 <- c(rep(-u[z], 1))
    }
    else if (lshape == "triangular") {
      if (ushape == "triangular") {
        l1 <- c(-u[z] * ((1 -3* ratio/max(ratio))/sqrt(ratio)))[1]
      }
    }
    
    if (ushape == "obf") {
      u2 <- (u[z] * 1/sqrt(ratio))[2]
    }
    else if (ushape == "pocock") {
      u2 <- rep(u[z], 1)
    }
    else if (ushape == "triangular") {
      u2 <- (u[z] * (1 + ratio/max(ratio))/sqrt(ratio))[2]
    }
    
    
    prob <-  1-pnorm(u1)+
      pmvnorm(lower = c(-Inf,u1,-Inf, u2,-Inf, -Inf),
              upper = c(u1,Inf,Inf, Inf,Inf,Inf),
              sigma = cov)[1]+
      pmvnorm(lower = c(l1,l1,-Inf, u2,-Inf, -Inf),
              upper = c(u1,u1,Inf,Inf,Inf,Inf),
              sigma = cov)[1]+
      pmvnorm(lower = c(l1,-Inf,-Inf, u2,-Inf, -Inf),
              upper = c(u1,l1,Inf,Inf, Inf, Inf),
              sigma = cov)[1]+
      pmvnorm(lower = c(-Inf,-Inf, u1,u2,-Inf,-Inf),
              upper = c(l1,u1,Inf,Inf,Inf,Inf),
              sigma = cov)[1]
    
    
    if(prob  <= alpha){
      alow <- u[z]
      finalalpha <- prob
      upperbounds <- c(u1,u2)
      lowerbounds <- l1
      
    }
    
    if(prob > alpha){
      
      aup <- u[z]
      break
      
    }
    
  }
  
  
  results <- list(alow,aup,upperbounds, lowerbounds,finalalpha)
  
  names(results) <- c("alow", "aup","upperbounds", "lowerbounds","totalpha")
  return(results)
}


### Function to compute SS in ORD NI design for 4 arms and 2 stages #####

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
#' @seed: seed for numerical integration
#' @power: type of power - "Reject all hyp." or "Reject at least one hyp." or "Reject LongANDMedium" 

boundaries_4arm2stageORDNonInf_Zmod <- function(theta,
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
                                                seed,
                                                power
){ 
  
  library(mvtnorm)
  library(cubature)
  
  pc <- theta[1]
  
  pt1 <- pc-noninf
  
  pt2 <- pc-noninf
  
  pt3 <- pc-noninf
  
  varpc <- pc*(1-pc)
  
  vart1 <- pt1*(1-pt1)
  
  vart2 <- pt2*(1-pt2)
  
  vart3 <- pt3*(1-pt3)
  
  var <- c(vart1, vart2, vart3)
  
  pnull <- varpc+var
  
  
  # 4-arm covariance matrix
  
  covbfin12 <- (sqrt(rhonull[1]*rhonull[2])*pc*(1-pc))/(sqrt(rhonull[4]*vart1+rhonull[1]*varpc)*sqrt(rhonull[4]*vart2+rhonull[2]*varpc))
  
  covbfin13 <- (sqrt(rhonull[1]*rhonull[3])*pc*(1-pc))/(sqrt(rhonull[4]*vart1+rhonull[1]*varpc)*sqrt(rhonull[4]*vart3+rhonull[3]*varpc))
  
  covbfin14 <- sqrt(rhonull[1]*rhonull[4]*r[2]*rhonull[1]*r[2]*rhonull[4])/sqrt((r[2]*rhonull[4]*vart1+r[2]*rhonull[1]*varpc)*(r[1]*rhonull[4]*vart1+r[1]*rhonull[1]*varpc))*((vart1/(r[2]*rhonull[1]))+(varpc/(r[2]*rhonull[4])))
  
  covbfin15 <- sqrt(rhonull[1]*rhonull[4]*r[2]*rhonull[2]/(r[2]*rhonull[4]))*varpc/(sqrt(rhonull[4]*vart1+rhonull[1]*varpc)*sqrt(r[2]*rhonull[4]*vart2+r[2]*rhonull[2]*varpc))
  
  covbfin16 <- sqrt(rhonull[1]*rhonull[4]*r[2]*rhonull[3]/(r[2]*rhonull[4]))*varpc/(sqrt(rhonull[4]*vart1+rhonull[1]*varpc)*sqrt(r[2]*rhonull[4]*vart3+r[2]*rhonull[3]*varpc))
  
  covbfin23 <- (sqrt(rhonull[2]*rhonull[3])*pc*(1-pc))/(sqrt(rhonull[4]*vart2+rhonull[2]*varpc)*sqrt(rhonull[4]*vart3+rhonull[3]*varpc))
  
  covbfin24 <- sqrt(rhonull[2]*rhonull[4]*r[2]*rhonull[1]/(r[2]*rhonull[4]))*varpc/(sqrt(rhonull[4]*vart2+rhonull[2]*varpc)*sqrt(r[2]*rhonull[4]*vart1+r[2]*rhonull[1]*varpc))
  
  covbfin25 <- sqrt(rhonull[2]*rhonull[4]*r[2]*rhonull[2]*r[2]*rhonull[4])/sqrt((r[2]*rhonull[4]*vart2+r[2]*rhonull[2]*varpc)*(r[1]*rhonull[4]*vart2+r[1]*rhonull[2]*varpc))*((vart2/(r[2]*rhonull[2]))+(varpc/(r[2]*rhonull[4])))
  
  covbfin26 <- sqrt(rhonull[2]*rhonull[4]*r[2]*rhonull[3]/(r[2]*rhonull[4]))*varpc/(sqrt(rhonull[4]*vart2+rhonull[2]*varpc)*sqrt(r[2]*rhonull[4]*vart3+r[2]*rhonull[3]*varpc))
  
  covbfin34 <- sqrt(rhonull[3]*rhonull[4]*r[2]*rhonull[1]/(r[2]*rhonull[4]))*varpc/(sqrt(rhonull[4]*vart3+rhonull[3]*varpc)*sqrt(r[2]*rhonull[4]*vart1+r[2]*rhonull[1]*varpc))
  
  covbfin35 <- sqrt(rhonull[3]*rhonull[4]*r[2]*rhonull[2]/(r[2]*rhonull[4]))*varpc/(sqrt(rhonull[4]*vart3+rhonull[3]*varpc)*sqrt(r[2]*rhonull[4]*vart2+r[2]*rhonull[2]*varpc))
  
  covbfin36 <- sqrt(rhonull[3]*rhonull[4]*r[2]*rhonull[3]*r[2]*rhonull[4])/sqrt((r[2]*rhonull[4]*vart3+r[2]*rhonull[3]*varpc)*(r[1]*rhonull[4]*vart3+r[1]*rhonull[3]*varpc))*((vart3/(r[2]*rhonull[3]))+(varpc/(r[2]*rhonull[4])))
  
  covbfin45 <- (sqrt(r[2]*rhonull[1]*r[2]*rhonull[2])*pc*(1-pc))/(sqrt(r[2]*rhonull[4]*vart1+r[2]*rhonull[1]*varpc)*sqrt(r[2]*rhonull[4]*vart2+r[2]*rhonull[2]*varpc))
  
  covbfin46 <- (sqrt(r[2]*rhonull[1]*r[2]*rhonull[3])*pc*(1-pc))/(sqrt(r[2]*rhonull[4]*vart1+r[2]*rhonull[1]*varpc)*sqrt(r[2]*rhonull[4]*vart3+r[2]*rhonull[3]*varpc))
  
  covbfin56 <-  (sqrt(r[2]*rhonull[2]*r[2]*rhonull[3])*pc*(1-pc))/(sqrt(r[2]*rhonull[4]*vart2+r[2]*rhonull[2]*varpc)*sqrt(r[2]*rhonull[4]*vart3+r[2]*rhonull[3]*varpc))
  
  covbfin <- matrix(c(1,covbfin12,covbfin13,covbfin14,covbfin15,covbfin16,
                      covbfin12, 1, covbfin23,covbfin24,covbfin25,covbfin26,
                      covbfin13,covbfin23,1,covbfin34,covbfin35,covbfin36,
                      covbfin14,covbfin24,covbfin34,1, covbfin45,covbfin46,
                      covbfin15,covbfin25,covbfin35, covbfin45,1, covbfin56,
                      covbfin16,covbfin26,covbfin36, covbfin46,covbfin56,1), nrow=6, byrow=TRUE)
  
  
  firststagematrix <- covbfin[1:(arms-1), 1:(arms-1)]
  
  
  u11 <- seq(from = 7, to = 0, by = -1)
  
  first <- bounds_4arm2stageORD(u11, arms,r, 
                                alpha, covbfin, firststagematrix,
                                ushape, lshape,prec)
  
  low <- first$alow
  
  up <- first$aup
  
  p <- rep(1, times = prec)
  
  for (p in 1:length(p)){
    
    callf <- bounds_4arm2stageORD(seq(from = low, 
                                      to = up, by = -1/(10^p)),
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
  
  
  totalalpha <- callf$totalpha
  
  u1 <- callf$upperbounds[1]
  
  u2 <- callf$upperbounds[2]
  
  l1 <- callf$lowerbounds[1]
  
  pop <- seq(1,4000, by = 1)
  
  maxsample <- NULL
  
  betat <- NULL
  
  ################ Power for 4arm-2stage design
  
  pc <- theta[1]
  
  pt1 <- theta[2]
  
  pt2 <- theta[3]
  
  pt3 <- theta[4]
  
  varpc <- pc*(1-pc)
  
  vart1 <- pt1*(1-pt1)
  
  vart2 <- pt2*(1-pt2)
  
  vart3 <- pt3*(1-pt3)
  
  var <- c(vart1, vart2, vart3)
  
  covbfin12 <- (sqrt(rhopower[1]*rhopower[2])*pc*(1-pc))/(sqrt(rhopower[4]*vart1+rhopower[1]*varpc)*sqrt(rhopower[4]*vart2+rhopower[2]*varpc))
  
  covbfin13 <- (sqrt(rhopower[1]*rhopower[3])*pc*(1-pc))/(sqrt(rhopower[4]*vart1+rhopower[1]*varpc)*sqrt(rhopower[4]*vart3+rhopower[3]*varpc))
  
  covbfin14 <- sqrt(rhopower[1]*rhopower[4]*r[2]*rhopower[1]*r[2]*rhopower[4])/sqrt((r[2]*rhopower[4]*vart1+r[2]*rhopower[1]*varpc)*(r[1]*rhopower[4]*vart1+r[1]*rhopower[1]*varpc))*((vart1/(r[2]*rhopower[1]))+(varpc/(r[2]*rhopower[4])))
  
  covbfin15 <- sqrt(rhopower[1]*rhopower[4]*r[2]*rhopower[2]/(r[2]*rhopower[4]))*varpc/(sqrt(rhopower[4]*vart1+rhopower[1]*varpc)*sqrt(r[2]*rhopower[4]*vart2+r[2]*rhopower[2]*varpc))
  
  covbfin16 <- sqrt(rhopower[1]*rhopower[4]*r[2]*rhopower[3]/(r[2]*rhopower[4]))*varpc/(sqrt(rhopower[4]*vart1+rhopower[1]*varpc)*sqrt(r[2]*rhopower[4]*vart3+r[2]*rhopower[3]*varpc))
  
  covbfin23 <- (sqrt(rhopower[2]*rhopower[3])*pc*(1-pc))/(sqrt(rhopower[4]*vart2+rhopower[2]*varpc)*sqrt(rhopower[4]*vart3+rhopower[3]*varpc))
  
  covbfin24 <- sqrt(rhopower[2]*rhopower[4]*r[2]*rhopower[1]/(r[2]*rhopower[4]))*varpc/(sqrt(rhopower[4]*vart2+rhopower[2]*varpc)*sqrt(r[2]*rhopower[4]*vart1+r[2]*rhopower[1]*varpc))
  
  covbfin25 <- sqrt(rhopower[2]*rhopower[4]*r[2]*rhopower[2]*r[2]*rhopower[4])/sqrt((r[2]*rhopower[4]*vart2+r[2]*rhopower[2]*varpc)*(r[1]*rhopower[4]*vart2+r[1]*rhopower[2]*varpc))*((vart2/(r[2]*rhopower[2]))+(varpc/(r[2]*rhopower[4])))
  
  covbfin26 <- sqrt(rhopower[2]*rhopower[4]*r[2]*rhopower[3]/(r[2]*rhopower[4]))*varpc/(sqrt(rhopower[4]*vart2+rhopower[2]*varpc)*sqrt(r[2]*rhopower[4]*vart3+r[2]*rhopower[3]*varpc))
  
  covbfin34 <- sqrt(rhopower[3]*rhopower[4]*r[2]*rhopower[1]/(r[2]*rhopower[4]))*varpc/(sqrt(rhopower[4]*vart3+rhopower[3]*varpc)*sqrt(r[2]*rhopower[4]*vart1+r[2]*rhopower[1]*varpc))
  
  covbfin35 <- sqrt(rhopower[3]*rhopower[4]*r[2]*rhopower[2]/(r[2]*rhopower[4]))*varpc/(sqrt(rhopower[4]*vart3+rhopower[3]*varpc)*sqrt(r[2]*rhopower[4]*vart2+r[2]*rhopower[2]*varpc))
  
  covbfin36 <- sqrt(rhopower[3]*rhopower[4]*r[2]*rhopower[3]*r[2]*rhopower[4])/sqrt((r[2]*rhopower[4]*vart3+r[2]*rhopower[3]*varpc)*(r[1]*rhopower[4]*vart3+r[1]*rhopower[3]*varpc))*((vart3/(r[2]*rhopower[3]))+(varpc/(r[2]*rhopower[4])))
  
  covbfin45 <- (sqrt(r[2]*rhopower[1]*r[2]*rhopower[2])*pc*(1-pc))/(sqrt(r[2]*rhopower[4]*vart1+r[2]*rhopower[1]*varpc)*sqrt(r[2]*rhopower[4]*vart2+r[2]*rhopower[2]*varpc))
  
  covbfin46 <- (sqrt(r[2]*rhopower[1]*r[2]*rhopower[3])*pc*(1-pc))/(sqrt(r[2]*rhopower[4]*vart1+r[2]*rhopower[1]*varpc)*sqrt(r[2]*rhopower[4]*vart3+r[2]*rhopower[3]*varpc))
  
  covbfin56 <-  (sqrt(r[2]*rhopower[2]*r[2]*rhopower[3])*pc*(1-pc))/(sqrt(r[2]*rhopower[4]*vart2+r[2]*rhopower[2]*varpc)*sqrt(r[2]*rhopower[4]*vart3+r[2]*rhopower[3]*varpc))
  
  covbfinB <- matrix(c(1,covbfin12,covbfin13,covbfin14,covbfin15,covbfin16,
                       covbfin12, 1, covbfin23,covbfin24,covbfin25,covbfin26,
                       covbfin13,covbfin23,1,covbfin34,covbfin35,covbfin36,
                       covbfin14,covbfin24,covbfin34,1, covbfin45,covbfin46,
                       covbfin15,covbfin25,covbfin35, covbfin45,1, covbfin56,
                       covbfin16,covbfin26,covbfin36, covbfin46,covbfin56,1), nrow=6, byrow=TRUE)
  
  for (n in 1:length(pop)){
    
    prob <- 0
    
    sigma1 <- (varpc/(rhopower[4]*pop[n]))+(vart1/(rhopower[1]*pop[n]))
    
    sigma2 <- (varpc/(rhopower[4]*pop[n]))+(vart2/(rhopower[2]*pop[n]))
    
    sigma3 <- (varpc/(rhopower[4]*pop[n]))+(vart3/(rhopower[3]*pop[n]))
    
    sigman <- c(sigma1,sigma2,sigma3)*rep(1/r, each=arms-1)
    
    mean1 <- rep(((theta[2:(arms)]-theta[1])+noninf),2)/(sqrt(sigman))
    
    if(power == "Reject all hyp."){
      
      set.seed(seed)
      
      prob <-pmvnorm(lower = c(rep(u1,3)),
                     upper = c(rep(Inf, 3)),
                     sigma = covbfinB[1:3,1:3],
                     mean = mean1[1:3])[1]+
        pmvnorm(lower = c(u1,u1,l1, -Inf, -Inf, u2),
                upper = c(Inf, Inf, u1, Inf, Inf, Inf),
                sigma = covbfinB,
                mean = mean1)[1]+
        pmvnorm(lower = c(u1,l1,l1, -Inf, u2, u2),
                upper = c(Inf, u1, u1, Inf, Inf, Inf),
                sigma = covbfinB,
                mean = mean1)[1]+
        pmvnorm(lower = c(u1,-Inf,u1, -Inf, u2, u2),
                upper = c(Inf, u1, Inf, Inf, Inf, Inf),
                sigma = covbfinB,
                mean = mean1)[1]+
        pmvnorm(lower = c(l1,l1,l1, u2, u2, u2),
                upper = c(u1, u1, u1, Inf, Inf, Inf),
                sigma = covbfinB,
                mean = mean1)[1]+
        pmvnorm(lower = c(-Inf,-Inf,u1, u2, u2, u2),
                upper = c(u1, Inf, Inf, Inf, Inf, Inf),
                sigma = covbfinB,
                mean = mean1)[1]+
        pmvnorm(lower = c(-Inf,u1,l1, u2, u2, u2),
                upper = c(u1, Inf, u1, Inf, Inf, Inf),
                sigma = covbfinB,
                mean = mean1)[1]
      
    }
    if(power == "Reject at least one hyp."){
      
      set.seed(seed)
      
      prob <- 1-pnorm(u1, sd = covbfinB[1,1],
                      mean = mean1[1])+
        pmvnorm(lower = c(-Inf,u1,-Inf, u2,-Inf, -Inf),
                upper = c(u1,Inf,Inf, Inf,Inf,Inf),
                sigma = covbfinB,
                mean = mean1)[1]+
        pmvnorm(lower = c(l1,l1,-Inf, u2,-Inf, -Inf),
                upper = c(u1,u1,Inf,Inf,Inf,Inf),
                sigma = covbfinB,
                mean = mean1)[1]+
        pmvnorm(lower = c(l1,-Inf,-Inf, u2,-Inf, -Inf),
                upper = c(u1,l1,Inf,Inf, Inf, Inf),
                sigma = covbfinB,
                mean = mean1)[1]+
        pmvnorm(lower = c(-Inf,-Inf, u1,u2,-Inf,-Inf),
                upper = c(l1,u1,Inf,Inf,Inf,Inf),
                sigma = covbfinB,
                mean = mean1)[1]
    }
    
    if(power == "Reject LongANDMedium"){
      
      set.seed(seed)
      
      prob <- pmvnorm(lower = c(rep(u1,2)),
                      upper = c(rep(Inf, 2)),
                      sigma = covbfinB[1:2,1:2],
                      mean = mean1[1:2])[1]+
        pmvnorm(lower = c(u1,l1,-Inf, -Inf, u2, -Inf),
                upper = c(Inf, u1, Inf, Inf, Inf, Inf),
                sigma = covbfinB,
                mean = mean1)[1]+
        pmvnorm(lower = c(u1,-Inf,u1, -Inf, u2, -Inf),
                upper = c(Inf, l1, Inf, Inf, Inf, Inf),
                sigma = covbfinB,
                mean = mean1)[1]+
        pmvnorm(lower = c(l1,l1,-Inf, u2, u2, -Inf),
                upper = c(u1, u1,Inf, Inf, Inf, Inf),
                sigma = covbfinB,
                mean = mean1)[1]+
        pmvnorm(lower = c(-Inf,u1,-Inf, u2, u2, -Inf),
                upper = c(u1, Inf, Inf, Inf, Inf, Inf),
                sigma = covbfinB,
                mean = mean1)[1]+
        pmvnorm(lower = c(-Inf,-Inf, u1,u2,u2,-Inf),
                upper = c(l1,u1,Inf,Inf,Inf,Inf),
                sigma = covbfinB,
                mean = mean1)[1]
    }
    
    
    if(prob > 1-beta){
      
      maxsample <- pop[n]
      
      totalbeta <- prob
      
      break
      
    }
    
    
  }
  
  results <- list(totalalpha, totalbeta,
                  c(u1,u2),l1,maxsample,callf$alow)
  
  names(results) <- c("sum_alphat", "totalbeta",
                      "upperbounds1", "lowerbounds1", "sample size per arm per stage",
                      "alow")
  
  return(results)
}


#### SIMULATION 4 arms - 2 stages ORD NI design ######

# PARAMETERS for both designs

seed <- 64736 # set seed to simulate the data

nsim <- 10^6 #number of simulations

arms <- 4 #number of arms

stage <- 2 # number of stages

alpha <- 0.05 # alpha-level of the test

ratio <- c(1,2) #allocation ratio sample size at different stages: n1=r[1], n2=r[2]*n1,...

noninf <- 0.1 # non-inferiority margin

pc <- 0.92 # response rate on the control arm

pt1 <- pc # response rate on the first treatment arm

pt2 <- pc # response rate on the second treatment arm

pt3 <- pc # response rate on the third treatment arm

beta <- 0.2 # beta-level 

rhonull <- rhopower <- c(1,1,1,1) #sample size in each treatment divided by sample size in the control: rho[1]=n11/n01, rho[2]=n21/n01, rho[3]=1

ushape <- "triangular" # shape of upper critical boundary

lshape <- "triangular" # shape of lower critical boundary

prec <- 3 # precision - number of digits after decimal point - for the critical boundaries

powerreq <- "Reject all hyp."#"Reject at least one hyp."#, "Reject LongANDMedium", "Reject all hyp."

# call function to find the sample size per arm per stage and the critical boundaries

fun <- boundaries_4arm2stageORDNonInf_Zmod(theta = c(pc,pt1,pt2,pt3),
                                      noninf = noninf,
                                      stage = stage,
                                      alpha = alpha,
                                      beta = beta,
                                      r = ratio,
                                      rhonull = rhonull,
                                      rhopower = rhopower,
                                      prec = prec,
                                      arms = arms,
                                      ushape = ushape,
                                      lshape = lshape,
                                      seed = seed,
                                      power=powerreq)
fun


### Scenario

####### RESTRUCTURE scenarios

plac <- rep(pc, 2)

trt1 <- c(pc-noninf, pc)

trt2 <- c(pc-noninf, pc)

trt3 <-  c(pc-noninf, pc)

scen <- cbind(plac,trt1,trt2, trt3)

scenario <- matrix(data = scen,
                   nrow = nrow(scen),
                   ncol = ncol(scen))

scenario <- cbind(scenario, scenario[,2]-scenario[,1], scenario[,3]-scenario[,1], scenario[,4]-scenario[,1])

colnames(scenario) <- c("p0", "p1", "p2", "p3","theta1", "theta2", "theta3")

scenario <- as.data.frame(scenario)

######## Sample size and bound grid values 

upper <- fun$alow#seq(round(fun$alow,3)-0.1,round(fun$alow,3)+0.1,by = 0.001)

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
    
    # Save rejections of rejecting H01 and H02 and H03
    
    rej_all <- list()
    
    rej_all[[1]] <- matrix(data=FALSE, nrow=nrow(scenario),
                           ncol=nsim)
    
    rej_all[[2]] <- matrix(data=FALSE, nrow=nrow(scenario),
                           ncol=nsim)
    
    # Save when rejecting any hypotheses
    
    rej_any <- list()
    
    rej_any[[1]] <- matrix(data=FALSE, nrow=nrow(scenario),
                           ncol=nsim)
    
    rej_any[[2]] <- matrix(data=FALSE, nrow=nrow(scenario),
                           ncol=nsim)
    
    # Save rejections of rejecting H01 and not H02,H03
    
    rej_ord_h01noh02noh03 <- matrix(data=FALSE, nrow=nrow(scenario),
                                    ncol=nsim)
    
    # Save rejections of rejecting H01 and H02 and not H03
    
    rej_ord_h01andh02noh03 <- matrix(data=FALSE, nrow=nrow(scenario),
                                     ncol=nsim)
    
    # Save all rejections for each scenario
    
    prop_ord <- list()
    
    prop_ord[[1]] <- matrix(data=NA, ncol=2,
                            nrow = nrow(scenario))
    prop_ord[[2]] <- matrix(data=NA, ncol=4,
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
  
  mean3 <- scenario[k,4]
  
  for (j in 1:nsim){
    
    plac1st <-  rbinom(n = popord, size = 1, prob = mean0)
    
    dur1st <-  rbinom(n = popord, size = 1, prob = mean1)
    
    dur2st <-  rbinom(n = popord, size = 1, prob = mean2)
    
    dur3st <-  rbinom(n = popord, size = 1, prob = mean3)
    
    estimatedss[k,j] <- length(plac1st)+length(dur1st)+length(dur2st)+length(dur3st)
    
    # Estimated response rates on the observed data
    
    p_plac11 <- mean(plac1st)
    
    p_dur11 <- mean(dur1st)
    
    p_dur12 <- mean(dur2st)
    
    p_dur13 <- mean(dur3st)
    
    var1 <- (mean1*(1-mean1)/length(dur1st))+(mean0*(1-mean0)/length(plac1st))
    
    var2 <- (mean2*(1-mean2)/length(dur2st))+(mean0*(1-mean0)/length(plac1st))
    
    var3 <- (mean3*(1-mean3)/length(dur3st))+(mean0*(1-mean0)/length(plac1st))
    
    # Z-statistics for the two treatment durations
    
    z1 <- (p_dur11-p_plac11+noninf)/sqrt(var1)
    
    z2 <- (p_dur12-p_plac11+noninf)/sqrt(var2)
    
    z3 <- (p_dur13-p_plac11+noninf)/sqrt(var3)
    
    # The null hypothesis is early rejected in the 1st stage
    
    if( z1 >= fun$upperbounds1[1] &  z2 >= fun$upperbounds1[1] 
        & z3 >= fun$upperbounds1[1]){
      
      rej_all[[1]][k,j] = TRUE
      
    }
    
    if(z1 <= fun$lowerbounds1[1] &  z2 <= fun$upperbounds1[1]  
       & z3 <= fun$upperbounds1[1]){
      
      rej_any[[1]][k,j] = TRUE
      
      
    }
    
    if(z1 >= fun$upperbounds1[1] &  z2 <= fun$lowerbounds1[1]  
       & z3 <= fun$upperbounds1[1]){
      
      rej_ord_h01noh02noh03[k,j] = TRUE
      
      
    }
    
    if(z1 >= fun$upperbounds1[1] &  z2 >= fun$upperbounds1[1]  
       & z3 <= fun$lowerbounds1[1]){
      
      rej_ord_h01andh02noh03[k,j] = TRUE
      
      
    }
    
    
    # Continue to second stage with all arms
    
    if((z1 <= fun$upperbounds1[1] &  z3 >=fun$upperbounds1[1])||
       (z1 <= fun$upperbounds1[1] &  z2 >=fun$upperbounds1[1] & 
        z3 >=fun$lowerbounds1[1] & z3 < fun$upperbounds1[1])||
       ( z1 > fun$lowerbounds1[1]  &  z1 < fun$upperbounds1[1] &
         z2 > fun$lowerbounds1[1]  &  z2 < fun$upperbounds1[1]
         & z3 > fun$lowerbounds1[1] & z3 < fun$upperbounds1[1] )
    ){
      
      # Simulate placebo and treatment populations
      
      plac21_2 <- rbinom(n = popord, size = 1, prob = mean0)
      
      dur21_2 <- rbinom(n = popord, size = 1, prob = mean1)
      
      dur22_2 <- rbinom(n = popord, size = 1, prob = mean2)
      
      dur23_2 <- rbinom(n = popord, size = 1, prob = mean3)
      
      estimatedss[k,j] <-  estimatedss[k,j]+length(plac21_2)+length(dur21_2)+length(dur22_2)+length(dur23_2)
      
      plac21 <- c(plac21_2, plac1st)
      
      dur21 <- c(dur21_2, dur1st)
      
      dur22 <- c(dur22_2, dur2st)
      
      dur23 <- c(dur23_2, dur3st)
      
      # Estimated response rates
      
      p_dur21 <- mean(dur21)
      
      p_plac21 <- mean(plac21)
      
      p_dur22 <- mean(dur22)
      
      p_dur23 <- mean(dur23)
      
      var1 <- (mean1*(1-mean1)/length(dur21))+(mean0*(1-mean0)/length(plac21))
      
      var2 <- (mean2*(1-mean2)/length(dur22))+(mean0*(1-mean0)/length(plac21))
      
      var3 <- (mean3*(1-mean3)/length(dur23))+(mean0*(1-mean0)/length(plac21))
      
      # Z-statistics for the two treatment durations
      
      z21 <- (p_dur21-p_plac21+noninf)/sqrt(var1)
      
      z22 <- (p_dur22-p_plac21+noninf)/sqrt(var2)
      
      z23 <- (p_dur23-p_plac21+noninf)/sqrt(var3)
      
      # Rejection of the null hypothesis at the second stage
      
      
      if( z21 >= fun$upperbounds1[2] &  
          z22 >= fun$upperbounds1[2] &
          z23 >= fun$upperbounds1[2]){
        
        rej_all[[2]][k,j] = TRUE
        
      }
      
      if( z21 < fun$upperbounds1[2]){
        
        rej_any[[2]][k,j] = TRUE
        
        
      }
      
      if( z21 >= fun$upperbounds1[2] &  
          z22 >= fun$upperbounds1[2] &
          z23 < fun$upperbounds1[2]){
        
        rej_ord_h01andh02noh03[k,j] = TRUE
        
      }
      
      if( z21 >= fun$upperbounds1[2] &  
          z22 < fun$upperbounds1[2]){
        
        rej_ord_h01noh02noh03[k,j] = TRUE
        
      }
      
    }
    # Continue to second stage with arms 1 and 2
    
    if((z1 <= fun$upperbounds1[1] &  z2 >=fun$upperbounds1[1] & z3 < fun$lowerbounds1[1] )||
       ( z1 > fun$lowerbounds1[1]  &  z1 < fun$upperbounds1[1] &
         z2 > fun$lowerbounds1[1]  &  z2 < fun$upperbounds1[1]
         & z3 < fun$lowerbounds1[1])
    ){
      
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
      
      var1 <- (mean1*(1-mean1)/length(dur21))+(mean0*(1-mean0)/length(plac21))
      
      var2 <- (mean2*(1-mean2)/length(dur22))+(mean0*(1-mean0)/length(plac21))
      
      # Z-statistics for the two treatment durations
      
      z21 <- (p_dur21-p_plac21+noninf)/sqrt(var1)
      
      z22 <- (p_dur22-p_plac21+noninf)/sqrt(var2)
      
      # Rejection of the null hypothesis at the second stage
      
      if (z21 < fun$upperbounds1[2]){
        
        rej_any[[2]][k,j]  = TRUE
        
      }
      
      if( z21 >= fun$upperbounds1[2] &  
          z22 >= fun$upperbounds1[2]){
        
        rej_ord_h01andh02noh03[k,j] = TRUE
        
      }
      
      if( z21 >= fun$upperbounds1[2] &  
          z22 < fun$upperbounds1[2] ){
        
        rej_ord_h01noh02noh03[k,j] = TRUE
        
      }
      
    }
    
    # Continue to second stage with arms 2 and 3
    
    if((z1 >=  fun$upperbounds1[1] &
        z2 > fun$lowerbounds1[1]  &  z2 < fun$upperbounds1[1]
        & z3 > fun$lowerbounds1[1]  &  z3 < fun$upperbounds1[1])||
       (z1 >=  fun$upperbounds1[1] &
        z2 < fun$upperbounds1[1]
        &  z3 >= fun$upperbounds1[1])
    ){
      
      # Simulate placebo and treatment populations
      
      plac21_2 <- rbinom(n = popord, size = 1, prob = mean0)
      
      dur22_2 <- rbinom(n = popord, size = 1, prob = mean2)
      
      dur23_2 <- rbinom(n = popord, size = 1, prob = mean3)
      
      estimatedss[k,j] <-  estimatedss[k,j]+length(plac21_2)+length(dur22_2)+length(dur23_2)
      
      plac21 <- c(plac21_2, plac1st)
      
      dur22 <- c(dur22_2, dur2st)
      
      dur23 <- c(dur23_2, dur3st)
      
      # Estimated response rates
      
      p_plac21 <- mean(plac21)
      
      p_dur22 <- mean(dur22)
      
      p_dur23 <- mean(dur23)
      
      var2 <- (mean2*(1-mean2)/length(dur22))+(mean0*(1-mean0)/length(plac21))
      
      var3 <- (mean3*(1-mean3)/length(dur23))+(mean0*(1-mean0)/length(plac21))
      
      # Z-statistics for the two treatment durations
      
      
      z22 <- (p_dur22-p_plac21+noninf)/sqrt(var2)
      
      z23 <- (p_dur23-p_plac21+noninf)/sqrt(var3)
      
      
      # Rejection of the null hypothesis at the second stage
      
      if (z23 >= fun$upperbounds1[2] & z22 >= fun$upperbounds1[2]){
        
        rej_all[[2]][k,j] = TRUE
        
        
      }
      
      if( z22 >= fun$upperbounds1[2] &
          z23 < fun$upperbounds1[2]){
        
        rej_ord_h01andh02noh03[k,j] = TRUE
        
      }
      
      if( z22 < fun$upperbounds1[2]){
        
        rej_ord_h01noh02noh03[k,j] = TRUE
        
      }
      
      
    }
    # Continue to second stage with 3
    
    if((z1 >=  fun$upperbounds1[1] &
        z2 >=  fun$upperbounds1[1] & 
        z3 > fun$lowerbounds1[1] &  z3 < fun$upperbounds1[1])
    ){
      
      # Simulate placebo and treatment populations
      
      plac21_2 <- rbinom(n = popord, size = 1, prob = mean0)
 
      dur23_2 <- rbinom(n = popord, size = 1, prob = mean3)
      
      estimatedss[k,j] <-  estimatedss[k,j]+length(plac21_2)+length(dur23_2)
      
      plac21 <- c(plac21_2, plac1st)
 
      dur23 <- c(dur23_2, dur3st)
      
      # Estimated response rates
     
      p_plac21 <- mean(plac21)
     
      p_dur23 <- mean(dur23)
   
      var3 <- (mean3*(1-mean3)/length(dur23))+(mean0*(1-mean0)/length(plac21))
      
      # Z-statistics for the two treatment durations
      
      z23 <- (p_dur23-p_plac21+noninf)/sqrt(var3)
      
      
      # Rejection of the null hypothesis at the second stage
      
      if (z23 >= fun$upperbounds1[2]){
        
        rej_all[[2]][k,j] = TRUE
        
      }
      
      if(z23 < fun$upperbounds1[2]){
        
        rej_ord_h01andh02noh03[k,j] = TRUE
        
      }
      
      
    }
    
    # Continue to second stage with 1
    
    if((z1 > fun$lowerbounds1[1]  &  z1 < fun$upperbounds1[1] &
        z2 <  fun$lowerbounds1[1] &  z3 < fun$upperbounds1[1])
    ){
      
      plac21_2 <- rbinom(n = popord, size = 1, prob = mean0)
      
      dur21_2 <- rbinom(n = popord, size = 1, prob = mean1)
     
      estimatedss[k,j] <-  estimatedss[k,j]+length(plac21_2)+length(dur21_2)
      
      plac21 <- c(plac21_2, plac1st)
      
      dur21 <- c(dur21_2, dur1st)
   
      # Estimated response rates
      
      p_dur21 <- mean(dur21)
      
      p_plac21 <- mean(plac21)
     
      var1 <- (mean1*(1-mean1)/length(dur21))+(mean0*(1-mean0)/length(plac21))

      # Z-statistics for the two treatment durations
      
      z21 <- (p_dur21-p_plac21+noninf)/sqrt(var1)

      
      # Rejection of the null hypothesis at the second stage
      
      if (z21 < fun$upperbounds1[2]){
        
        rej_any[[2]][k,j] = TRUE
        
      }
      
      if( z21 >= fun$upperbounds1[2]){
        
        rej_ord_h01noh02noh03[k,j] = TRUE
        
      }
      
    }
    
    # Continue to second stage with 2
    
    if((z1 >= fun$upperbounds1[1] &
        z2 <  fun$upperbounds1[1] &  z2 > fun$lowerbounds1[1]
        &  z3 < fun$lowerbounds1[1])
    ){
      
      plac21_2 <- rbinom(n = popord, size = 1, prob = mean0)
      
      dur22_2 <- rbinom(n = popord, size = 1, prob = mean2)
   
      estimatedss[k,j] <-  estimatedss[k,j]+length(plac21_2)+length(dur22_2)
      
      plac21 <- c(plac21_2, plac1st)
     
      dur22 <- c(dur22_2, dur2st)
     
      # Estimated response rates

      p_plac21 <- mean(plac21)
      
      p_dur22 <- mean(dur22)
    
      var2 <- (mean2*(1-mean2)/length(dur22))+(mean0*(1-mean0)/length(plac21))
 
      # Z-statistics for the two treatment durations
      
      z22 <- (p_dur22-p_plac21+noninf)/sqrt(var2)
 
      # Rejection of the null hypothesis at the second stage
      
      
      if( z22 >= fun$upperbounds1[2]){
        
        rej_ord_h01andh02noh03[k,j] = TRUE
        
      }
      
      if( z22 < fun$upperbounds1[2]){
        
        rej_ord_h01noh02noh03[k,j] = TRUE
        
      }
      
    }
    
    
  }
  # Average probabilities of rej the hypotheses
  
  prop_ord[[1]][k,1] <- sum(rej_all[[1]][k,])/nsim
  
  prop_ord[[1]][k,2] <- sum(rej_any[[1]][k,])/nsim
  
  prop_ord[[2]][k,1] <- sum(rej_all[[2]][k,])/nsim
  
  prop_ord[[2]][k,2] <- sum(rej_any[[2]][k,])/nsim
  
  prop_ord[[2]][k,3] <- sum(rej_ord_h01andh02noh03[k,])/nsim
  
  prop_ord[[2]][k,4] <- sum(rej_ord_h01noh02noh03[k,])/nsim
  
  estimatedsamplesize[k] <- mean(estimatedss[k,])
  
}

    colnames(prop_ord[[1]]) <- c("ORD_rejall_1stage",
                                 "ORD_rejany_1stage")
    
    colnames(prop_ord[[2]]) <- c("ORD_rejall_2stage",
                                 "ORD_rejany_2stage",
                                 "ORD_rejH01andH02noH03",
                                 "ORD_rejH01noH02noH03")
    
    colnames(estimatedsamplesize) <- "ESS Ordered ORD"
    
    summary_ORD <- cbind(popord,
                         ceiling(popord*2*arms),
                         scenario,
                         prop_ord,
                         estimatedsamplesize,
                         u1,
                         u2,
                         l1)
    
    colnames(summary_ORD)[1:2]<- c("patients_perarm1stage ORD",
                                   "total sample ORD")
    
    summary_ORD <- summary_ORD %>% mutate(
      
      rej_all_overall = ORD_rejall_1stage + ORD_rejall_2stage,
      
      rej_any_overall = ORD_rejany_1stage + ORD_rejany_2stage,
      
      rejLM = ORD_rejH01andH02noH03+ORD_rejall_1stage + ORD_rejall_2stage,
      
      typeI = ORD_rejH01noH02noH03 + ORD_rejH01andH02noH03 + rej_all_overall,
      
      sumprob = typeI + rej_any_overall
      
    )
    
    typeI <- summary_ORD[1,]$`Reject at least one hyp.`
    
    power <- summary_ORD[2,powerreq]
    
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




# 4-arm 2-stage MAMS(m) design

rm(list=ls())

library(dplyr)

##### Time to first positive claim

# Durations of the treatment

durations <- c(6,4,3,2)

# Recruitment rate per month

recrate <- 30

#### Function to find critical bounds in 4-arm 2-stage MAMS(m) NI design #####

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
#' @lfix: fixed value for lower critical bound

bounds_4arm2stagem <- function(u1, interimN, arms,r,
                               alpha, cov, first, ushape,lshape,lfix, prec){
  
  library(mvtnorm)
  library(gtools)
  
  alphat <- rep(0,length = interimN)
  
  alphatcaseup<- NULL
  
  alphatcasedown<- NULL
  
  alphatcasefin<- NULL
  
  for(j in 1:length(u1)){
    
    for (i in 1:interimN){
      
      if (ushape == "obf") {
        u <- u1[j] * 1/sqrt(r)
      }
      else if (ushape == "pocock") {
        u <- rep(u1[j], i)
      }
      else if (ushape == "triangular") {
        u <- u1[j] * (1 + r/max(r))/sqrt(r)
      }
      
      if (lshape == "obf") {
        l <- c(-u1[j] * 1/sqrt(r))[1:(i-1)]
      }
      else if (lshape == "pocock") {
        l <- c(rep(-u1[j], i-1 ))
      }
      else if (lshape == "triangular") {
        if (ushape == "triangular") {
          l <- c(-u1[j] * ((1 -3* r/max(r))/sqrt(r)))[1:(i-1)]
        }
      }
      else if (lshape == "fixed") {
        l <- c(rep(lfix, i-1 ))
      }
      
      
      if(i ==1){
        set.seed(123456)
        
        alphat[i] <- 1-pmvnorm(lower = rep(-Inf,times = arms-1),
                               upper = rep(u[1],times = arms-1),
                               sigma = first)[1]
        
      }
      if(i>1){
        
        set.seed(123456)
        
        alphat[i] <- pmvnorm(lower = c(l[1],-Inf, -Inf, u[2], -Inf, -Inf),
                             upper = c(u[1],l[1], l[1], Inf, Inf, Inf),
                             sigma = cov)[1]+
          pmvnorm(lower = c(-Inf, l[1],-Inf,  -Inf,u[2], -Inf),
                  upper = c(l[1],u[1], l[1], Inf, Inf, Inf),
                  sigma = cov)[1]+
          pmvnorm(lower = c(-Inf, -Inf,l[1],  -Inf, -Inf,u[2]),
                  upper = c(l[1],l[1], u[1], Inf, Inf, Inf),
                  sigma = cov)[1]+
          pmvnorm(lower = c(l[1],l[1], -Inf, u[2], -Inf, -Inf),
                  upper = c(u[1],u[1], l[1], Inf, Inf, Inf),
                  sigma = cov)[1]+
          pmvnorm(lower = c(l[1],l[1], -Inf, -Inf, u[2], -Inf),
                  upper = c(u[1],u[1], l[1], u[2], Inf, Inf),
                  sigma = cov)[1]+
          pmvnorm(lower = c(l[1], -Inf,l[1], u[2], -Inf, -Inf),
                  upper = c(u[1], l[1],u[1], Inf, Inf, Inf),
                  sigma = cov)[1]+
          pmvnorm(lower = c(l[1], -Inf,l[1], -Inf, -Inf,u[2]),
                  upper = c(u[1],l[1], u[1], u[2], Inf, Inf),
                  sigma = cov)[1]+
          pmvnorm(lower = c( -Inf,l[1],l[1],  -Inf,u[2], -Inf),
                  upper = c(l[1], u[1],u[1], Inf, Inf, Inf),
                  sigma = cov)[1]+
          pmvnorm(lower = c( -Inf,l[1],l[1],  -Inf,-Inf,u[2]),
                  upper = c(l[1],u[1], u[1],  Inf,u[2], Inf),
                  sigma = cov)[1]+
          pmvnorm(lower = c(l[1],l[1], l[1], u[2], -Inf, -Inf),
                  upper = c(u[1],u[1], u[1], Inf, Inf, Inf),
                  sigma = cov)[1]+
          pmvnorm(lower = c(l[1],l[1], l[1],  -Inf,u[2], -Inf),
                  upper = c(u[1],u[1], u[1], Inf, Inf, Inf),
                  sigma = cov)[1]+
          pmvnorm(lower = c(l[1],l[1], l[1],  -Inf, -Inf,u[2]),
                  upper = c(u[1],u[1], u[1], Inf, Inf, Inf),
                  sigma = cov)[1]-
          pmvnorm(lower = c(l[1],l[1], l[1], u[2], u[2], -Inf),
                  upper = c(u[1],u[1], u[1], Inf, Inf, Inf),
                  sigma = cov)[1]-
          pmvnorm(lower = c(l[1],l[1], l[1], u[2],  -Inf,u[2]),
                  upper = c(u[1],u[1], u[1], Inf, Inf, Inf),
                  sigma = cov)[1]-
          pmvnorm(lower = c(l[1],l[1], l[1],  -Inf,u[2],u[2]),
                  upper = c(u[1],u[1], u[1], Inf, Inf, Inf),
                  sigma = cov)[1]+
          pmvnorm(lower = c(l[1],l[1], l[1], u[2],u[2],u[2]),
                  upper = c(u[1],u[1], u[1], Inf, Inf, Inf),
                  sigma = cov)[1]
      }
      
      
    }
    
    if(sum(alphat) <= alpha){
      
      if(i==1){
        
        upperbound <- u[1]
        lowerbound <- u[1]
        totalpha <- sum(alphat)
        alow <- u1[j]
        finalalpha <- alphat
        
      }
      else{
        upperbound <- u
        lowerbound <- l
        totalpha <- sum(alphat)
        alow <- u1[j]
        finalalpha <- alphat
      }
      
    }
    
    if(sum(alphat) > alpha){
      
      aup <- u1[j]
      break
      
      
    }
    
  }
  
  results <- list(upperbound,alow, aup, lowerbound,totalpha)
  
  names(results) <- c("upperbound", "alow", "aup", "lowerbound", "totalpha")
  return(results)
  
}


##### Function to compute ss in a 4-arm 2-stage MAMS(m) NI design #####

#' @theta: vector of clinically relevant difference 
#' @noninf: non-inferiority margin 
#' @interimN: number of interim analyses
#' @alpha: alpha-level of the hypothesis test
#' @beta: beta-level of the test
#' @r: allocation ratio sample size in the first and second stage
#' @rhonull: allocation ratio vector for the null configuration - (r_1^(1), r_1^(2), ..., r_1^(0))
#' @rhopower: allocation ratio vector for the power configuration - (r_1^(1), r_1^(2), ..., r_1^(0))
#' @prec: precision for the boundaries: numbers after the comma
#' @arms: number of arms (3)
#' @ushape: "pocock", "obf", "triangular"
#' @lshape: "pocock", "obf", "triangular"
#' @lfix: fixed value for lower critical bound
#' @power: type of power - "reject all" or "reject at least one" or "reject L&M" 

boundaries_4arm2stage <- function(theta,
                                  noninf, 
                                  interimN, 
                                  alpha, 
                                  beta,
                                  r, 
                                  rhonull,
                                  rhopower, 
                                  prec,
                                  arms,
                                  ushape,
                                  lshape,
                                  lfix,
                                  power){ 
  
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
  
  
  
  u1 <- seq(from = 7, to = 0, by = -1)
  
  first <- bounds_4arm2stagem(u1, stage, arms,r, 
                              alpha, covbfin, firststagematrix,
                              ushape, lshape,lfix,prec)
  
  low <- first$alow
  
  up <- first$aup
  
  p <- rep(1, times = prec)
  
  for (p in 1:length(p)){
    
    callf <- bounds_4arm2stagem(seq(from = low, 
                                    to = up, by = -1/(10^p)), 
                                interimN= stage,
                                arms,
                                r,
                                alpha,
                                covbfin,
                                firststagematrix,
                                ushape,
                                lshape,
                                lfix,
                                prec)
    
    low <- callf$alow
    
    up <- callf$aup
  }
  
  
  upperbound <- callf$upperbound
  lowerbound <- callf$lowerbound
  a <- callf$alow
  totalalpha <- callf$totalpha
  
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
  
  firststagematrixB <- covbfinB[1:(arms-1), 1:(arms-1)]
  
  
  maxsample <- NULL
  
  betat <- NULL
  
  betatcaseup <- NULL
  
  betatcasedown <- NULL
  
  betatcasefin <- NULL
  
  pop <- seq(from = 1, to = 500, by = 1)
  
  
  for (n in 1:length(pop)){
    
    prob <- 0
    
    sigma1 <- (varpc/(rhopower[4]*pop[n]))+(vart1/(rhopower[1]*pop[n]))
    
    sigma2 <- (varpc/(rhopower[4]*pop[n]))+(vart2/(rhopower[2]*pop[n]))
    
    sigma3 <- (varpc/(rhopower[4]*pop[n]))+(vart3/(rhopower[3]*pop[n]))
    
    sigman <- c(sigma1,sigma2,sigma3)*rep(1/r, each=arms-1)
    
    mean1 <- rep(((theta[2:(arms)]-theta[1])+noninf),2)/(sqrt(sigman))
    
    for (i in 1:interimN){
      
      if(power == "reject all"){
        
        #hard code for 4-arm 2-stage MAMS(m) design
        
        if(arms==4){
          
          if(i==1){
            betat[i] <- pmvnorm(lower = rep(upperbound[1],times = arms-1),
                                upper = rep(Inf,times = arms-1),
                                sigma = firststagematrixB,
                                mean = mean1[1:3])[1]
          }
          
          if(i==2){
            l <- lowerbound
            u <- upperbound
            
            betat[i] <-  pmvnorm(lower = c(l[1], l[1],l[1], u[2], u[2], u[2]),
                                 upper = c(u[1],u[1],u[1], Inf, Inf, Inf),
                                 sigma = covbfinB,
                                 mean = mean1)[1]+
              pmvnorm(lower = c(l[1], u[1],u[1], u[2], -Inf,-Inf),
                      upper = c(u[1],Inf,Inf, Inf, Inf, Inf),
                      sigma = covbfinB,
                      mean = mean1)[1]+
              pmvnorm(lower = c(u[1], l[1],u[1], -Inf, u[2], -Inf),
                      upper = c(Inf,u[1],Inf, Inf, Inf, Inf),
                      sigma = covbfinB,
                      mean = mean1)[1]+
              pmvnorm(lower = c(u[1], u[1],l[1],-Inf, -Inf, u[2]),
                      upper = c(Inf,Inf,u[1], Inf, Inf, Inf),
                      sigma = covbfinB,
                      mean = mean1)[1]+
              pmvnorm(lower = c(l[1], l[1],u[1], u[2], u[2], -Inf),
                      upper = c(u[1],u[1],Inf, Inf, Inf, Inf),
                      sigma = covbfinB,
                      mean = mean1)[1]+
              pmvnorm(lower = c(u[1], l[1],l[1], -Inf, u[2], u[2]),
                      upper = c(Inf,u[1],u[1], Inf, Inf, Inf),
                      sigma = covbfinB,
                      mean = mean1)[1]+
              pmvnorm(lower = c(l[1], u[1],l[1], u[2],-Inf, u[2]),
                      upper = c(u[1],Inf,u[1], Inf, Inf, Inf),
                      sigma = covbfinB,
                      mean = mean1)[1]
          }
        }
        
        
      }
      
      
      if(power == "reject at least one"){
        
        if(i ==1){
          set.seed(123456)
          
          betat[i] <- 1-pmvnorm(lower = rep(-Inf,times = arms-1),
                                upper = rep(upperbound[1],times = arms-1),
                                sigma = firststagematrixB,
                                mean = mean1[1:3])[1]
          
        }
        if(i>1){
          
          l <- lowerbound
          u <- upperbound
          
          betat[i] <- pmvnorm(lower = c(l[1],-Inf, -Inf, u[2], -Inf, -Inf),
                              upper = c(u[1],l[1], l[1], Inf, Inf, Inf),
                              sigma = covbfinB,
                              mean = mean1)[1]+
            pmvnorm(lower = c(-Inf, l[1],-Inf,  -Inf,u[2], -Inf),
                    upper = c(l[1],u[1], l[1], Inf, Inf, Inf),
                    sigma = covbfinB,
                    mean = mean1)[1]+
            pmvnorm(lower = c(-Inf, -Inf,l[1],  -Inf, -Inf,u[2]),
                    upper = c(l[1],l[1], u[1], Inf, Inf, Inf),
                    sigma = covbfinB,
                    mean = mean1)[1]+
            pmvnorm(lower = c(l[1],l[1], -Inf, u[2], -Inf, -Inf),
                    upper = c(u[1],u[1], l[1], Inf, Inf, Inf),
                    sigma = covbfinB,
                    mean = mean1)[1]+
            pmvnorm(lower = c(l[1],l[1], -Inf, -Inf, u[2], -Inf),
                    upper = c(u[1],u[1], l[1], u[2], Inf, Inf),
                    sigma = covbfinB,
                    mean = mean1)[1]+
            pmvnorm(lower = c(l[1], -Inf,l[1], u[2],-Inf,-Inf),
                    upper = c(u[1], l[1],u[1], Inf, Inf, Inf),
                    sigma = covbfinB,
                    mean = mean1)[1]+
            pmvnorm(lower = c(l[1], -Inf,l[1], -Inf, -Inf,u[2]),
                    upper = c(u[1],l[1], u[1], u[2], Inf, Inf),
                    sigma = covbfinB,
                    mean = mean1)[1]+
            pmvnorm(lower = c( -Inf,l[1],l[1],  -Inf,u[2], -Inf),
                    upper = c(l[1], u[1],u[1], Inf, Inf, Inf),
                    sigma = covbfinB,
                    mean = mean1)[1]+
            pmvnorm(lower = c( -Inf,l[1],l[1],  -Inf,-Inf,u[2]),
                    upper = c(l[1],u[1], u[1],  Inf,u[2], Inf),
                    sigma = covbfinB,
                    mean = mean1)[1]+
            pmvnorm(lower = c(l[1],l[1], l[1], u[2], -Inf, -Inf),
                    upper = c(u[1],u[1], u[1], Inf, Inf, Inf),
                    sigma = covbfinB,
                    mean = mean1)[1]+
            pmvnorm(lower = c(l[1],l[1], l[1],  -Inf,u[2], -Inf),
                    upper = c(u[1],u[1], u[1], Inf, Inf, Inf),
                    sigma = covbfinB,
                    mean = mean1)[1]+
            pmvnorm(lower = c(l[1],l[1], l[1],  -Inf, -Inf,u[2]),
                    upper = c(u[1],u[1], u[1], Inf, Inf, Inf),
                    sigma = covbfinB,
                    mean = mean1)[1]-
            pmvnorm(lower = c(l[1],l[1], l[1], u[2], u[2], -Inf),
                    upper = c(u[1],u[1], u[1], Inf, Inf, Inf),
                    sigma = covbfinB,
                    mean = mean1)[1]-
            pmvnorm(lower = c(l[1],l[1], l[1], u[2],  -Inf,u[2]),
                    upper = c(u[1],u[1], u[1], Inf, Inf, Inf),
                    sigma = covbfinB,
                    mean = mean1)[1]-
            pmvnorm(lower = c(l[1],l[1], l[1],  -Inf,u[2],u[2]),
                    upper = c(u[1],u[1], u[1], Inf, Inf, Inf),
                    sigma = covbfinB,
                    mean = mean1)[1]+
            pmvnorm(lower = c(l[1],l[1], l[1], u[2],u[2],u[2]),
                    upper = c(u[1],u[1], u[1], Inf, Inf, Inf),
                    sigma = covbfinB,
                    mean = mean1)[1]
          
          
        }
      }
      
      if(power == "reject L&M"){
        
        if(i ==1){
          set.seed(123456)
          
          betat[i] <- pmvnorm(lower = c(upperbound[1],upperbound[1], -Inf),
                              upper = c(Inf,Inf, Inf),
                              sigma = firststagematrixB,
                              mean = mean1[1:3])[1]
          
        }
        if(i>1){
          
          l <- lowerbound
          u <- upperbound
          
          betat[i] <- pmvnorm(lower = c(l[1],l[1], -Inf, u[2], u[2], -Inf),
                              upper = c(u[1],u[1], Inf, Inf, Inf, Inf),
                              sigma = covbfinB,
                              mean = mean1)[1]+
            pmvnorm(lower = c(u[1],l[1], -Inf, -Inf, u[2], -Inf),
                    upper = c(Inf,u[1], Inf, Inf, Inf, Inf),
                    sigma = covbfinB,
                    mean = mean1)[1]+
            pmvnorm(lower = c(l[1],u[1],-Inf, u[2], -Inf, -Inf),
                    upper = c(u[1],Inf, Inf, Inf, Inf, Inf),
                    sigma = covbfinB,
                    mean = mean1)[1]
          
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
                      "lowerbounds", "upperbounds", "sample size per arm per stage", "a")
  
  return(results)
}


#### SIMULATION 4-arm 2-stage MAMS(m) NI design ######

# PARAMETERS for both designs

seed <- 64736 # set seed to simulate the data

nsim <- 5*10^4 #number of simulations

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

rho <-rhonull <- rhopower <- c(1,1,1,1) #sample size in each treatment divided by sample size in the control: rho[1]=n11/n01, rho[2]=n21/n01, rho[3]=1

ushape <- "triangular" # shape of upper critical boundary

lshape <- "triangular" # shape of lower critical boundary

prec <- 3 # precision - number of digits after decimal point - for the critical boundaries

powereq <- "reject all" # "reject at least one", "reject L&M"

# call function to find the sample size per arm per stage and the critical boundaries

interim <- 0.5 #at which proportion of the first sample size we do the interim analysis --> at 40% of the tot population 

func <- boundaries_4arm2stage(theta = c(pc,pt1,pt2,pt3),
                              noninf = noninf,
                              interimN = stage,
                              alpha = alpha,
                              beta = beta,
                              r = ratio,
                              rhonull=rho,
                              rhopower=rho,
                              prec = prec,
                              arms = arms,
                              ushape = ushape,
                              lshape =ushape,
                              power = powereq)

func 

####### RESTRUCTURE scenarios

# probability to allocate every patient on each treatment arm

p <- rep(1/arms, arms)

### Scenario

### RESTRUCTURE scenarios

plac <- c(0,	0,	0,	0,	0,	0,	0,	0,	0)+pc

trt1 <- c(0,	0,	0,	0,	0,	0,	0,	0,	-0.1)+pc

trt2 <- c(-0.1,	-0.07,	-0.04,	-0.01,	0.02,	0.05,	0.08,	0,	-0.1)+pc

trt3 <-  c(-0.11,	-0.08,	-0.06,	-0.03,	0,	0.02,	0.05,	0,	-0.1)+pc



scen <- cbind(plac,trt1,trt2, trt3)

scenario <- matrix(data = scen,
                   nrow = nrow(scen),
                   ncol = ncol(scen))

scenario <- cbind(scenario, scenario[,2]-scenario[,1], scenario[,3]-scenario[,1], scenario[,4]-scenario[,1])

colnames(scenario) <- c("p0", "p1", "p2", "p3","theta1", "theta2", "theta3")

scenario <- as.data.frame(scenario)

######## Sample size and bound grid values 

popord <-104


# triangular bounds

# critical boundaries

func$upperbounds[1] <- 2.361
func$lowerbounds[1] <- 0.787
func$upperbounds[2] <- 2.226

summary_MAMS <- NULL


# Save rejections of H01 and H02 not 3

rej_h012not3 <- list()

rej_h012not3[[1]] <- matrix(data=FALSE, nrow=nrow(scenario),
                            ncol=nsim)

rej_h012not3[[2]] <- matrix(data=FALSE, nrow=nrow(scenario),
                            ncol=nsim)

# Save rejections of H01 and H03 not 2

rej_h013not2 <- list()

rej_h013not2[[1]] <- matrix(data=FALSE, nrow=nrow(scenario),
                            ncol=nsim)

rej_h013not2[[2]] <- matrix(data=FALSE, nrow=nrow(scenario),
                            ncol=nsim)

# Save rejections of H02 and H03 not 1

rej_h023not1 <- list()

rej_h023not1[[1]] <- matrix(data=FALSE, nrow=nrow(scenario),
                            ncol=nsim)

rej_h023not1[[2]] <- matrix(data=FALSE, nrow=nrow(scenario),
                            ncol=nsim)

# Save rejections of H01 and H02 and H03

rej_h01andh02andh03 <- list()

rej_h01andh02andh03[[1]] <- matrix(data=FALSE, nrow=nrow(scenario),
                                   ncol=nsim)

rej_h01andh02andh03[[2]] <- matrix(data=FALSE, nrow=nrow(scenario),
                                   ncol=nsim)

# Save rejections of H01 not 2and3

rej_h01not23 <- list()

rej_h01not23[[1]] <- matrix(data=FALSE, nrow=nrow(scenario),
                            ncol=nsim)

rej_h01not23[[2]] <- matrix(data=FALSE, nrow=nrow(scenario),
                            ncol=nsim)

# Save rejections of H02 not 1and3

rej_h02not13 <- list()

rej_h02not13[[1]] <- matrix(data=FALSE, nrow=nrow(scenario),
                            ncol=nsim)

rej_h02not13[[2]] <- matrix(data=FALSE, nrow=nrow(scenario),
                            ncol=nsim)

# Save rejections of H03 not 1and2

rej_h03not12 <- list()

rej_h03not12[[1]] <- matrix(data=FALSE, nrow=nrow(scenario),
                            ncol=nsim)

rej_h03not12[[2]] <- matrix(data=FALSE, nrow=nrow(scenario),
                            ncol=nsim)


bothineff_all  <- matrix(data=FALSE, nrow=nrow(scenario),
                         ncol=nsim)

prop_rej <- list()

prop_rej[[1]] <- matrix(data=NA, ncol=7,
                        nrow = nrow(scenario))
prop_rej[[2]] <- matrix(data=NA, ncol=8,
                        nrow = nrow(scenario))

# Estimated sample size for each scenario. --> The average of the sample size used for each simulation

estimatedss <- actualmaxss <- matrix(data=NA, nrow=nrow(scenario),
                                     ncol=nsim)

estimatedsamplesize <- matrix(data=NA, ncol=1,
                              nrow = nrow(scenario))

durationtrial <- matrix(data=NA, ncol=1,
                        nrow = nrow(scenario))

timetofirstposclaim <- matrix(data=NA, ncol=1,
                              nrow = nrow(scenario))

dur <- matrix(data=0, nrow=nrow(scenario),
              ncol=nsim)

durfin <- matrix(data=0, nrow=nrow(scenario),
                 ncol=nsim)

tfpc <- matrix(data=0, nrow=nrow(scenario),
               ncol=nsim)

for (k in 1:nrow(scenario)){
  
  set.seed(seed = seed)
  
  mean0 <- scenario[k,1]
  
  mean1 <- scenario[k,2]
  
  mean2 <- scenario[k,3]
  
  mean3 <- scenario[k,4]
  
  for (j in 1:nsim){
    
    plac1st <- dur1st <- dur2st <- dur3st <- NULL
    
    plac21 <- dur21 <- dur22 <- dur23 <- NULL
    
    plac_2missing <- plac <- short <- long <- medium <- 0
    
    count <- NULL
    
    lastpatplac <- FALSE # when I recruit the last patient on the placebo
    
    lastpatlong <- FALSE # when I recruit the last patient on the longest
    
    lastpatshort <- FALSE # when I recruit the last patient on the shortest
    
    lastpatmedium <- FALSE # when I recruit the last patient on the medium
    
    #recruit until get popord on the shortest
    
    while(plac < popord || short < popord || long < popord || medium < popord){
      count <- c(count,sample(c("C","L","M", "S"),recrate, prob = p, replace=TRUE))
      short <- sum(count=="S")
      medium <- sum(count=="M")
      long <- sum(count=="L")
      plac <- sum(count=="C")
      
      if(plac<popord){
        lastpatplac <- TRUE
        month_lastpatplac <- length(count)/recrate
      }
      if(short<popord){
        lastpatshort <- TRUE
        month_lastpatshort <- length(count)/recrate
      }
      if(long<popord){
        lastpatlong <- TRUE
        month_lastpatlong <- length(count)/recrate
      }
      if(medium<popord){
        lastpatmedium <- TRUE
        month_lastpatmedium <- length(count)/recrate
      }
    }
    #table(count)
    
    n_plac <- sum(count=="C")
    n_long <- sum(count=="L")
    n_short <- sum(count=="S")
    n_medium <- sum(count=="M")
    months <- length(count)/recrate
    
    n_interim <- min(c(n_plac, n_long,n_short, n_medium))
    
    plac11st <- rbinom(n = n_interim, prob= mean0, size = 1)
    
    dur11st <- rbinom(n = n_interim, prob= mean1, size = 1)
    
    dur21st <- rbinom(n = n_interim, prob= mean2, size = 1)
    
    dur31st <- rbinom(n = n_interim, prob= mean3, size = 1)
    
    
    # continue to recruit until observe the patients on the shortest duration
    
    plac_over <- n_plac
    
    count_over <- NULL
    
    lastpatientsobserved <- c(month_lastpatplac, month_lastpatlong,
                              month_lastpatmedium,month_lastpatshort)
    
    time_interim_month <- max(c(month_lastpatplac+durations[1],#+6,
                                month_lastpatlong+durations[2],#+6,
                                month_lastpatmedium+durations[3],#+6,
                                month_lastpatshort+durations[4]))#+6))
    
    index_timeinter <- which.max(c(month_lastpatplac+durations[1],#+6,
                                   month_lastpatlong+durations[2],#+6,
                                   month_lastpatmedium+durations[3],#+6,
                                   month_lastpatshort+durations[4]))#+6))
    
    continuerecr <- time_interim_month-lastpatientsobserved[index_timeinter]
    
    dur[k,j] <- time_interim_month
    
    #recruit until get popord on the shortest
    
    while((length(count_over)/recrate) < continuerecr & plac_over < 2*popord){
      count_over <- c(count_over,sample(c("C","L","M","S"), recrate, prob = p, replace=TRUE))
      months_over <- length(count_over)/recrate
      plac_over <- sum(count_over=="C")
    }
    #table(count_over)
    
    n_plac_over <- sum(count_over=="C")
    n_long_over <- sum(count_over=="L")
    n_short_over <- sum(count_over=="S")
    n_medium_over <- sum(count_over=="M")
    
    plac1st <- c(plac11st, rbinom(n = (n_plac-n_interim)+n_plac_over, prob= mean0, size = 1))
    
    dur1st <- c(dur11st, rbinom(n = (n_long-n_interim)+n_long_over, prob= mean1, size = 1))
    
    dur2st <- c(dur21st, rbinom(n = (n_medium-n_interim)+n_medium_over, prob= mean2, size = 1))
    
    dur3st <- c(dur31st, rbinom(n = (n_short-n_interim)+n_short_over, prob= mean3, size = 1))
    
    estimatedss[k,j] <- length(plac1st)+length(dur1st)+length(dur2st)+length(dur3st)
    
    # Estimated response rates on the observed data
    
    p_plac11 <- mean(plac11st)
    
    p_dur11 <- mean(dur11st)
    
    p_dur12 <- mean(dur21st)
    
    p_dur13 <- mean(dur31st)
    
    p1 <- (mean1*(1-mean1)/length(dur11st))+(mean0*(1-mean0)/length(plac11st))
    
    p2 <- (mean2*(1-mean2)/length(dur21st))+(mean0*(1-mean0)/length(plac11st))
    
    p3 <- (mean3*(1-mean3)/length(dur31st))+(mean0*(1-mean0)/length(plac11st))
    
    sdpooled11 <- p1
    
    sdpooled12 <- p2
    
    sdpooled13 <- p3
    
    # Z-statistics for the two treatment durations
    
    z1 <- (p_dur11-p_plac11+noninf)/sqrt(sdpooled11)
    
    z2 <- (p_dur12-p_plac11+noninf)/sqrt(sdpooled12)
    
    z3 <- (p_dur13-p_plac11+noninf)/sqrt(sdpooled13)
    
    # Count when null hypothesis are rejected at the end of the stage
    
    
    if(
      (z1 >= func$upperbounds[1] & z2 >= func$upperbounds[1] &  z3 >= func$upperbounds[1])){ 
      
      
      rej_h01andh02andh03[[1]][k,j] = TRUE
      tfpc[k,j] = time_interim_month
      
    }
    
    
    if((z1 >= func$upperbounds[1] & z3 < func$lowerbounds[1] &  z2 < func$lowerbounds[1])){
      
      
      rej_h01not23[[1]][k,j] = TRUE
      tfpc[k,j] = time_interim_month
      
    }
    
    if((z2 >= func$upperbounds[1] & z1 < func$lowerbounds[1] &  z3 < func$lowerbounds[1])){
      
      
      rej_h02not13[[1]][k,j] = TRUE
      tfpc[k,j] = time_interim_month
      
    }
    
    if((z3 >= func$upperbounds[1] & z1 < func$lowerbounds[1] &  z2 < func$lowerbounds[1])){
      
      
      rej_h03not12[[1]][k,j] = TRUE
      tfpc[k,j] = time_interim_month
      
    }
    
    if((z1 >= func$upperbounds[1] & z3 < func$lowerbounds[1] &  z2 >= func$upperbounds[1])){
      
      
      rej_h012not3[[1]][k,j] = TRUE
      tfpc[k,j] = time_interim_month
      
    }
    
    if((z2 >= func$upperbounds[1] & z1 < func$lowerbounds[1] &  z3 >= func$upperbounds[1])){
      
      
      rej_h023not1[[1]][k,j] = TRUE
      tfpc[k,j] = time_interim_month
      
    }
    
    if((z1 >= func$upperbounds[1] & z2 < func$lowerbounds[1] &  z3 >= func$upperbounds[1])){
      
      
      rej_h013not2[[1]][k,j] = TRUE
      tfpc[k,j] = time_interim_month
      
    }
    
    if (z1 < func$lowerbounds[1] & z2 <func$lowerbounds[1] & z3 < func$lowerbounds[1]){
      
      
      bothineff_all[k,j] = TRUE
      tfpc[k,j] = time_interim_month
      
      
      
    }
    
    # Continue to the second stage with all arms
    
    if(z1 < func$upperbounds[1] & z1 > func$lowerbounds[1] &
       z2 < func$upperbounds[1] & z2 > func$lowerbounds[1] &
       z3 < func$upperbounds[1] & z3 > func$lowerbounds[1] ){
      
      plac21 <- plac1st
      dur21 <- dur1st
      dur22 <- dur2st
      dur23 <- dur3st
      
      
      # continue to recruit until observe the patients on the control duration
      
      plac_2missing <- if_else((2*popord-(length(plac21)))>0, 
                               (2*popord-(length(plac21))), -1)
      
      plac_2 <- months_2 <- 0
      
      count_2 <- NULL
      
      #recruit until get popord on the shortest
      
      while(plac_2 < plac_2missing){
        count_2 <- c(count_2,sample(c("C","L","M", "S"), recrate, prob = p, replace=TRUE))
        months_2 <- length(count_2)/recrate
        plac_2 <- sum(count_2=="C")
      }
      #table(count_2)
      
      n_plac_2 <- sum(count_2=="C")
      n_long_2 <- sum(count_2=="L")
      n_short_2 <- sum(count_2=="S")
      n_medium_2 <- sum(count_2=="M")
      
      # Simulate placebo and treatment populations
      
      plac21 <- c(plac21, rbinom(n = n_plac_2, prob= mean0, size = 1))
      
      dur21 <- c(dur21, rbinom(n = n_long_2, prob= mean1, size = 1))
      
      dur22 <- c(dur22, rbinom(n = n_medium_2, prob= mean2, size = 1))
      
      dur23 <- c(dur23, rbinom(n = n_short_2, prob= mean3, size = 1))
      
      #duration
      
      dur[k,j] = dur[k,j]+(months_2 + durations[1])*(months_2>0)+(0)*(months_2==0)#+6
      
      #ESS
      
      estimatedss[k,j] <-  length(plac21)+length(dur21)+length(dur22)+length(dur23)
      
      actualmaxss[k,j] <-  length(plac21)+length(dur21)+length(dur22)+length(dur23)
      
      # Estimated response rates
      
      p_dur21 <- mean(dur21)
      
      p_plac21 <- mean(plac21)
      
      p_dur22 <- mean(dur22)
      
      p_dur23 <- mean(dur23)
      
      p21 <- (mean1*(1-mean1)/length(dur21))+(mean0*(1-mean0)/length(plac21))
      
      p22 <- (mean2*(1-mean2)/length(dur22))+(mean0*(1-mean0)/length(plac21))
      
      p23 <- (mean3*(1-mean3)/length(dur23))+(mean0*(1-mean0)/length(plac21))
      
      sdpooled21 <- p21
      
      sdpooled22 <- p22
      
      sdpooled23 <- p23
      
      # Z-statistics for the two treatment durations
      
      z21 <- (p_dur21-p_plac21+noninf)/sqrt(sdpooled21)
      
      z22 <- (p_dur22-p_plac21+noninf)/sqrt(sdpooled22)
      
      z23 <- (p_dur23-p_plac21+noninf)/sqrt(sdpooled23)
      
      # Rejection of the null hypothesis at the second stage
      
      
      
      
      if((z22 < func$upperbounds[2] & z23 < func$upperbounds[2] &  z21 >= func$upperbounds[2])) {
        
        rej_h01not23[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
      }
      
      if(
        (z22 >= func$upperbounds[2] & z23 < func$upperbounds[2] &  z21 >= func$upperbounds[2])){
        
        rej_h012not3[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
      } 
      if(
        (z22 < func$upperbounds[2] & z23 >= func$upperbounds[2] &  z21 >= func$upperbounds[2])){
        
        rej_h013not2[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
      }
      
      if(
        (z21 >= func$upperbounds[2] & z22 >= func$upperbounds[2] &  z23 >= func$upperbounds[2])
      ){ 
        
        
        rej_h01andh02andh03[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
        
      }
      
      
      if((z21 < func$upperbounds[2] & z23 < func$upperbounds[2] &  z22 >= func$upperbounds[2])){
        
        rej_h02not13[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
      }
      
      if( 
        (z21 < func$upperbounds[2] & z23 >= func$upperbounds[2] &  z22 >= func$upperbounds[2]) ){
        
        rej_h023not1[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
      }
      
      
      # rej H03
      
      
      if((z21 < func$upperbounds[2] & z22 < func$upperbounds[2] &  z23 >= func$upperbounds[2])){
        
        rej_h03not12[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
      }
      
      
      
      
      
      if (z21 < func$upperbounds[2] & z22 <func$upperbounds[2] & z23 < func$upperbounds[2]){
        
        
        bothineff_all[k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
        
        
      }
      
      
      
    }
    
    
    # # Continue to the second stage with first and second arm
    
    if((z1 < func$upperbounds[1] & z1 >=  func$lowerbounds[1] &
        z2 < func$upperbounds[1] & z2 >=  func$lowerbounds[1] &
        z3 >=  func$upperbounds[1]) ||
       (z1 < func$upperbounds[1] & z1 >=  func$lowerbounds[1] &
        z2 < func$upperbounds[1] & z2 >=  func$lowerbounds[1] &
        z3 <=  func$lowerbounds[1])){
      
      # Simulate placebo and treatment populations
      
      plac21 <- plac1st
      dur21 <- dur1st
      dur22 <- dur2st
      
      # continue to recruit until observe the patients on the control duration
      
      plac_2missing <-  if_else((2*popord-(length(plac21)))>0, 
                                (2*popord-(length(plac21))), -1)
      
      plac_2 <-  months_2 <-0
      
      count_2 <- NULL
      
      #recruit until get popord on the shortest
      
      while(plac_2 < plac_2missing){
        count_2 <- c(count_2,sample(c("C","L","M"), recrate, prob = p[1:3], replace=TRUE))
        months_2 <- length(count_2)/recrate
        plac_2 <- sum(count_2=="C")
      }
      #table(count_2)
      
      n_plac_2 <- sum(count_2=="C")
      n_long_2 <- sum(count_2=="L")
      n_medium_2 <- sum(count_2=="M")
      
      # Simulate placebo and treatment populations
      
      plac21 <- c(plac21, rbinom(n = n_plac_2, prob= mean0, size = 1))
      
      dur21 <- c(dur21, rbinom(n = n_long_2, prob= mean1, size = 1))
      
      dur22 <- c(dur22, rbinom(n = n_medium_2, prob= mean2, size = 1))
      
      
      #duration
      
      dur[k,j] = dur[k,j]+(months_2 + durations[1])*(months_2>0)+(0)*(months_2==0)#+6
      
      #ESS
      
      estimatedss[k,j] <-  length(plac21)+length(dur21)+length(dur22)+length(dur3st)
      
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
      
      if( ( z21 >= func$upperbounds[2] & z22 < func$upperbounds[2] & z3 >= func$upperbounds[1]) ){
        
        rej_h013not2[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
        
        
      }
      
      if( ( z21 >= func$upperbounds[2] & z22 < func$upperbounds[2] & z3 < func$lowerbounds[1]) ){
        
        rej_h01not23[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
        
        
      }
      
      if( ( z22 >= func$upperbounds[2] & z21 < func$upperbounds[2] & z3 >= func$upperbounds[1]) ){
        
        rej_h023not1[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
        
        
      }
      
      if( ( z22 >= func$upperbounds[2] & z21 < func$upperbounds[2] & z3< func$lowerbounds[1]) ){
        
        rej_h02not13[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
        
        
      }
      
      if(( z21 >= func$upperbounds[2] & z22 >= func$upperbounds[2]) & z3 >= func$upperbounds[1]){
        
        rej_h01andh02andh03[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
        
      }
      
      if(( z21 >= func$upperbounds[2] & z22 >= func$upperbounds[2]) & z3 < func$lowerbounds[1]){
        
        rej_h012not3[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
        
      }
      
      
      if( ( z21 < func$upperbounds[2] & z22 < func$upperbounds[2]) & z3 >= func$upperbounds[1]){
        
        rej_h03not12[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
        
        
      }
      
      if( ( z21 < func$upperbounds[2] & z22 < func$upperbounds[2]) & z3 < func$lowerbounds[1]){
        
        bothineff_all[k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
        
        
      }
      
      
      
      
      
    }
    
    
    # Continue to the second stage with first and third arm
    
    if((z1 < func$upperbounds[1] & z1 >=  func$lowerbounds[1] &
        z3 < func$upperbounds[1] & z3 >=  func$lowerbounds[1] &
        z2 >=  func$upperbounds[1]) ||
       (z1 < func$upperbounds[1] & z1 >=  func$lowerbounds[1] &
        z3 < func$upperbounds[1] & z3 >=  func$lowerbounds[1] &
        z2 <  func$lowerbounds[1])){
      
      plac21 <- plac1st
      dur21 <- dur1st
      
      dur23 <- dur3st
      
      # continue to recruit until observe the patients on the control duration
      
      plac_2missing <-  if_else((2*popord-(length(plac21)))>0, 
                                (2*popord-(length(plac21))), -1)
      
      plac_2 <- months_2 <-0
      
      count_2 <- NULL
      
      #recruit until get popord on the shortest
      
      while(plac_2 < plac_2missing){
        count_2 <- c(count_2,sample(c("C","L","S"),  recrate, prob = p[c(1,2,4)], replace=TRUE))
        months_2 <- length(count_2)/recrate
        plac_2 <- sum(count_2=="C")
      }
      #table(count_2)
      
      n_plac_2 <- sum(count_2=="C")
      n_long_2 <- sum(count_2=="L")
      n_short_2 <- sum(count_2=="S")
      
      
      # Simulate placebo and treatment populations
      
      plac21 <- c(plac21, rbinom(n = n_plac_2, prob= mean0, size = 1))
      
      dur21 <- c(dur21, rbinom(n = n_long_2, prob= mean1, size = 1))
      
      dur23 <- c(dur23, rbinom(n = n_short_2, prob= mean3, size = 1))
      
      #duration
      
      dur[k,j] = dur[k,j]+(months_2 + durations[1])*(months_2>0)+(0)*(months_2==0)#+6
      
      #ESS
      
      estimatedss[k,j] <-  length(plac21)+length(dur21)+length(dur23)+length(dur2st)
      
      # Estimated response rates
      
      p_dur21 <- mean(dur21)
      
      p_plac21 <- mean(plac21)
      
      
      p_dur23 <- mean(dur23)
      
      p21 <- (mean1*(1-mean1)/length(dur21))+(mean0*(1-mean0)/length(plac21))
      
      
      p23 <- (mean3*(1-mean3)/length(dur23))+(mean0*(1-mean0)/length(plac21))
      
      sdpooled21 <- p21
      
      
      sdpooled23 <- p23
      
      # Z-statistics for the two treatment durations
      
      z21 <- (p_dur21-p_plac21+noninf)/sqrt(sdpooled21)
      
      
      z23 <- (p_dur23-p_plac21+noninf)/sqrt(sdpooled23)
      
      
      if( ( z21 >= func$upperbounds[2] & z23 < func$upperbounds[2] & z2 >=  func$upperbounds[1]) ){
        
        rej_h012not3[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
        
        
      }
      
      if( ( z21 >= func$upperbounds[2] & z23 < func$upperbounds[2] & z2 <  func$lowerbounds[1]) ){
        
        rej_h01not23[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
        
        
      }
      
      if( ( z23 >= func$upperbounds[2] & z21 < func$upperbounds[2] & z2 >=  func$upperbounds[1]) ){
        
        rej_h023not1[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
        
        
      }
      
      if( ( z23 >= func$upperbounds[2] & z21 < func$upperbounds[2] & z2 <  func$lowerbounds[1]) ){
        
        rej_h03not12[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
        
        
      }
      
      if(
        ( z21 >= func$upperbounds[2] & z23 >= func$upperbounds[2]  & z2 >=  func$upperbounds[1])){
        
        
        rej_h01andh02andh03[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
        
      }
      
      if(
        ( z21 >= func$upperbounds[2] & z23 >= func$upperbounds[2]  & z2<  func$lowerbounds[1])){
        
        
        rej_h013not2[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
        
      }
      
      
      if( ( z21 < func$upperbounds[2] & z23 < func$upperbounds[2]& z2 >=  func$upperbounds[1]) ){
        
        rej_h02not13[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
        
        
      }
      
      if( ( z21 < func$upperbounds[2] & z23 < func$upperbounds[2]& z2 <  func$lowerbounds[1]) ){
        
        bothineff_all[k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
        
        
      }
      
      
    }
    
    # Continue to the second stage with third and second arm
    
    if((z3 < func$upperbounds[1] & z3 >=  func$lowerbounds[1] &
        z2 < func$upperbounds[1] & z2 >=  func$lowerbounds[1] &
        z1 >=  func$upperbounds[1])||
       (z3 < func$upperbounds[1] & z3 >=  func$lowerbounds[1] &
        z2 < func$upperbounds[1] & z2 >=  func$lowerbounds[1] &
        z1 <  func$lowerbounds[1])){
      
      # Simulate placebo and treatment populations
      
      plac21 <- plac1st
      
      dur22 <- dur2st
      dur23 <- dur3st
      
      # continue to recruit until observe the patients on the control duration
      
      plac_2missing <- if_else((2*popord-(length(plac21)))>0, 
                               (2*popord-(length(plac21))), -1)
      
      plac_2 <- months_2 <- 0
      
      count_2 <- NULL
      
      #recruit until get popord on the shortest
      
      while(plac_2 < plac_2missing){
        count_2 <- c(count_2,sample(c("C","M", "S"), recrate, prob = p[c(1,3,4)], replace=TRUE))
        months_2 <- length(count_2)/recrate
        plac_2 <- sum(count_2=="C")
      }
      #table(count_2)
      
      n_plac_2 <- sum(count_2=="C")
      
      n_short_2 <- sum(count_2=="S")
      n_medium_2 <- sum(count_2=="M")
      
      
      # Simulate placebo and treatment populations
      
      plac21 <- c(plac21, rbinom(n = n_plac_2, prob= mean0, size = 1))
      
      dur22 <- c(dur22, rbinom(n = n_medium_2, prob= mean2, size = 1))
      
      dur23 <- c(dur23, rbinom(n = n_short_2, prob= mean3, size = 1))
      
      #duration
      
      dur[k,j] =   dur[k,j]+(months_2 + durations[1])*(months_2>0)+(0)*(months_2==0)#+6
      
      #ESS
      
      estimatedss[k,j] <-  length(plac21)+length(dur22)+length(dur23)+length(dur1st)
      
      # Estimated response rates
      
      p_dur23 <- mean(dur23)
      
      p_plac21 <- mean(plac21)
      
      p_dur22 <- mean(dur22)
      
      p22 <- (mean2*(1-mean2)/length(dur22))+(mean0*(1-mean0)/length(plac21))
      
      p23 <- (mean3*(1-mean3)/length(dur23))+(mean0*(1-mean0)/length(plac21))
      
      
      sdpooled22 <- p22
      
      sdpooled23 <- p23
      
      # Z-statistics for the two treatment durations
      
      
      z22 <- (p_dur22-p_plac21+noninf)/sqrt(sdpooled22)
      
      z23 <- (p_dur23-p_plac21+noninf)/sqrt(sdpooled23)
      
      if( ( z22 >= func$upperbounds[2] & z23 < func$upperbounds[2] &
            z1 >=  func$upperbounds[1]) ){
        
        rej_h012not3[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
        
        
      }
      
      if( ( z22 >= func$upperbounds[2] & z23 < func$upperbounds[2] &
            z1 <  func$lowerbounds[1]) ){
        
        rej_h023not1[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
        
        
      }
      
      if( ( z23 >=func$upperbounds[2] & z22 < func$upperbounds[2] &
            z1 >=  func$upperbounds[1])  ){
        
        rej_h013not2[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
        
        
      }
      
      if( ( z23 >=func$upperbounds[2] & z22 < func$upperbounds[2] &
            z1 <  func$lowerbounds[1])){
        
        rej_h03not12[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
        
        
      }
      
      if( 
        ( z22 >= func$upperbounds[2] & z23 >= func$upperbounds[2]&
          z1 >=  func$upperbounds[1])){
        
        
        rej_h01andh02andh03[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
        
      }
      
      if( 
        ( z22 >= func$upperbounds[2] & z23 >= func$upperbounds[2]&
          z1 <  func$lowerbounds[1])){
        
        
        rej_h023not1[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
        
      }
      
      
      if( ( z22 < func$upperbounds[2] & z23 < func$upperbounds[2]&
            z1 >=  func$upperbounds[1]) ){
        
        rej_h01not23[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
        
        
      }
      
      if( ( z22 < func$upperbounds[2] & z23 < func$upperbounds[2]&
            z1<  func$lowerbounds[1]) ){
        
        bothineff_all[k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
        
        
      }
      
      
      
      
      
    }
    
    # Continue to the second stage with first arm
    
    if((z2 >= func$upperbounds[1] &
        z3 >=  func$upperbounds[1] & z1 < func$upperbounds[1] & z1 > func$lowerbounds[1]) ||
       (z2 < func$lowerbounds[1] &
        z3 < func$lowerbounds[1] & z1 < func$upperbounds[1] & z1 > func$lowerbounds[1]) ||
       (z2 >= func$upperbounds[1] &
        z3 < func$lowerbounds[1] & z1 < func$upperbounds[1] & z1 > func$lowerbounds[1]) ||
       (z2 < func$lowerbounds[1] &
        z3 >= func$upperbounds[1] & z1 < func$upperbounds[1] & z1 > func$lowerbounds[1])){
      
      plac21 <- plac1st
      dur21 <- dur1st
      
      # continue to recruit until observe the patients on the control duration
      
      plac_2missing <- if_else((2*popord-(length(plac21)))>0, 
                               (2*popord-(length(plac21))), -1)
      
      plac_2 <- months_2 <-0
      
      count_2 <- NULL
      
      #recruit until get popord on the shortest
      
      while(plac_2 < plac_2missing){
        count_2 <- c(count_2,sample(c("C","L"), recrate, prob = p[1:2], replace=TRUE))
        months_2 <- length(count_2)/recrate
        plac_2 <- sum(count_2=="C")
      }
      #table(count_2)
      
      n_plac_2 <- sum(count_2=="C")
      n_long_2 <- sum(count_2=="L")
      
      
      # Simulate placebo and treatment populations
      
      plac21 <- c(plac21, rbinom(n = n_plac_2, prob= mean0, size = 1))
      
      dur21 <- c(dur21, rbinom(n = n_long_2, prob= mean1, size = 1))
      
      #duration
      
      dur[k,j] = dur[k,j]+(months_2 + durations[1])*(months_2>0)+(0)*(months_2==0)#+6
      
      #ESS
      
      
      estimatedss[k,j] <- length(plac21)+length(dur21)+length(dur3st)+length(dur2st)
      
      
      # Estimated response rates
      
      p_plac21 <- mean(plac21)
      
      p_dur21 <- mean(dur21)
      
      p21 <- (mean1*(1-mean1)/length(dur21))+(mean0*(1-mean0)/length(plac21))
      
      sdpooled21 <- p21
      
      # Z-statistics for the two treatment durations
      
      z21 <- (p_dur21-p_plac21+noninf)/sqrt(sdpooled21)
      
      
      # rej H01
      
      if( z21 >= func$upperbounds[2] & z2 >= func$upperbounds[1] &
          z3 >=  func$upperbounds[1] ){
        
        
        rej_h01andh02andh03[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
      }
      
      if( z21 >= func$upperbounds[2] & z2 < func$lowerbounds[1] &
          z3 <  func$lowerbounds[1] ){
        
        
        rej_h01not23[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
      }
      
      if( z21 >= func$upperbounds[2] & z2 >= func$upperbounds[1] &
          z3 <  func$lowerbounds[1] ){
        
        
        rej_h012not3[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
      }
      
      if( z21 >= func$upperbounds[2] & z3 >= func$upperbounds[1] &
          z2 <  func$lowerbounds[1] ){
        
        
        rej_h013not2[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
      }
      
      if( z21 < func$upperbounds[2] & z2 >= func$upperbounds[1] &
          z3 >=  func$upperbounds[1]){
        
        
        rej_h023not1[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
      }
      
      if( z21 < func$upperbounds[2] & z2 >= func$upperbounds[1] &
          z3 <  func$lowerbounds[1]){
        
        
        rej_h02not13[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
      }
      
      if( z21 < func$upperbounds[2] & z3 >= func$upperbounds[1] &
          z2 <  func$lowerbounds[1]){
        
        
        rej_h03not12[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
      }
      
      if( z21 < func$upperbounds[2] & z2 < func$lowerbounds[1] &
          z3 <  func$lowerbounds[1] ){
        
        
        bothineff_all[k,j] = TRUE
        tfpc[k,j] = dur[k,j]
      }
      
      
      
    }
    
    # Continue to the second stage with second arm
    
    if((z1 >= func$upperbounds[1] &
        z3 >=  func$upperbounds[1] & z2 < func$upperbounds[1] & z2 >= func$lowerbounds[1]) ||
       (z1 < func$lowerbounds[1] &
        z3 < func$lowerbounds[1] & z2 < func$upperbounds[1] & z2 >= func$lowerbounds[1]) ||
       (z1 >= func$upperbounds[1] &
        z3 < func$lowerbounds[1] & z2 < func$upperbounds[1] & z2 > func$lowerbounds[1]) ||
       (z1 < func$lowerbounds[1] &
        z3 >= func$upperbounds[1] & z2 < func$upperbounds[1] & z2 > func$lowerbounds[1])){
      
      plac21 <- plac1st
      
      dur22 <- dur2st
      
      
      # continue to recruit until observe the patients on the control duration
      
      plac_2missing <-  if_else((2*popord-(length(plac21)))>0, 
                                (2*popord-(length(plac21))), -1)
      
      plac_2 <- months_2 <-0
      
      count_2 <- NULL
      
      #recruit until get popord on the shortest
      
      while(plac_2 < plac_2missing){
        count_2 <- c(count_2,sample(c("C","M"), recrate, prob = p[c(1,3)], replace=TRUE))
        months_2 <- length(count_2)/recrate
        plac_2 <- sum(count_2=="C")
      }
      #table(count_2)
      
      n_plac_2 <- sum(count_2=="C")
      
      n_medium_2 <- sum(count_2=="M")
      
      # Simulate placebo and treatment populations
      
      plac21 <- c(plac21, rbinom(n = n_plac_2, prob= mean0, size = 1))
      
      dur22 <- c(dur22, rbinom(n = n_medium_2, prob= mean2, size = 1))
      
      #duration
      
      dur[k,j] = dur[k,j]+(months_2 + durations[1])*(months_2>0)+(0)*(months_2==0)#+6
      
      #ESS
      
      estimatedss[k,j] <- length(plac21)+length(dur22)+length(dur1st)+length(dur3st)
      
      # Estimated response rates
      
      p_plac21 <- mean(plac21)
      
      p_dur22 <- mean(dur22)
      
      p22 <- (mean2*(1-mean2)/length(dur22))+(mean0*(1-mean0)/length(plac21))
      
      sdpooled22 <- p22
      
      # Z-statistics for the two treatment durations
      
      z22 <- (p_dur22-p_plac21+noninf)/sqrt(sdpooled22)
      
      # rej H02
      
      if( z22 >= func$upperbounds[2] & z1 >= func$upperbounds[1] &
          z3 >=  func$upperbounds[1]){
        
        rej_h01andh02andh03[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
      }
      
      if( z22 >= func$upperbounds[2] & z1 < func$lowerbounds[1] &
          z3 < func$lowerbounds[1]){
        
        rej_h02not13[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
      }
      
      if( z22 >= func$upperbounds[2] & z1 >= func$upperbounds[1] &
          z3 < func$lowerbounds[1]){
        
        rej_h012not3[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
      }
      
      if( z22 >= func$upperbounds[2] & z3 >= func$upperbounds[1] &
          z1 < func$lowerbounds[1]){
        
        rej_h023not1[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
      }
      
      if( z22 < func$upperbounds[2]& z1 >= func$upperbounds[1] &
          z3 >=  func$upperbounds[1]){
        
        rej_h013not2[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
      }
      
      if( z22 < func$upperbounds[2]& z1 < func$lowerbounds[1] &
          z3 >=  func$upperbounds[1]){
        
        rej_h03not12[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
      }
      
      if( z22 < func$upperbounds[2]& z3 < func$lowerbounds[1] &
          z1 >=  func$upperbounds[1]){
        
        rej_h01not23[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
      }
      
      if( z22 < func$upperbounds[2]& z1 < func$lowerbounds[1] &
          z3 < func$lowerbounds[1]){
        
        bothineff_all[k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
      }
      
      
    }
    
    # Continue to the second stage with third arm
    
    if((z1 >= func$upperbounds[1] &
        z2 >=  func$upperbounds[1] & z3 < func$upperbounds[1] & z3 >= func$lowerbounds[1]) ||
       (z1 < func$lowerbounds[1] &
        z2 < func$lowerbounds[1] & z3 < func$upperbounds[1] & z3 >= func$lowerbounds[1]) ||
       (z1 >= func$upperbounds[1] &
        z2 < func$lowerbounds[1] & z3 < func$upperbounds[1] & z3 > func$lowerbounds[1]) ||
       (z1 < func$lowerbounds[1] &
        z2 >= func$upperbounds[1] & z3 < func$upperbounds[1] & z3 > func$lowerbounds[1])){
      
      # Simulate placebo and treatment populations
      
      plac21 <- plac1st
      
      dur23 <- dur3st
      
      # continue to recruit until observe the patients on the control duration
      
      plac_2missing <- if_else((2*popord-(length(plac21)))>0, 
                              (2*popord-(length(plac21))), -1)
      
      plac_2 <- months_2 <-0
      
      count_2 <- NULL
      
      #recruit until get popord on the shortest
      
      while(plac_2 < plac_2missing){
        count_2 <- c(count_2,sample(c("C","S"), recrate, prob = p[c(1,4)], replace=TRUE))
        months_2 <- length(count_2)/recrate
        plac_2 <- sum(count_2=="C")
      }
      #table(count_2)
      
      n_plac_2 <- sum(count_2=="C")
      
      n_short_2 <- sum(count_2=="S")
      
      
      # Simulate placebo and treatment populations
      
      plac21 <- c(plac21, rbinom(n = n_plac_2, prob= mean0, size = 1))
      
      dur23 <- c(dur23, rbinom(n = n_short_2, prob= mean3, size = 1))
      
      #duration
      
      dur[k,j] = dur[k,j]+(months_2 + durations[1])*(months_2>0)+(0)*(months_2==0)#+6
      
      #ESS
      
      
      estimatedss[k,j] <- length(plac21)+length(dur23)+length(dur1st)+length(dur2st)
      
      
      # Estimated response rates
      
      p_plac21 <- mean(plac21)
      
      p_dur23 <- mean(dur23)
      
      
      p23 <- (mean3*(1-mean3)/length(dur23))+(mean0*(1-mean0)/length(plac21))
      
      
      sdpooled23 <- p23
      
      # Z-statistics for the two treatment durations
      
      z23 <- (p_dur23-p_plac21+noninf)/sqrt(sdpooled23)
      
      
      # rej H03
      
      if( z23 >= func$upperbounds[2] & z1 >= func$upperbounds[1] &
          z2 >=  func$upperbounds[1]){
        
        rej_h01andh02andh03[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
      }
      
      if( z23 >= func$upperbounds[2] & z1 < func$lowerbounds[1] &
          z2 < func$lowerbounds[1]){
        
        rej_h03not12[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
      }
      
      if( z23 >= func$upperbounds[2] & z1 >= func$upperbounds[1] &
          z2 < func$lowerbounds[1]){
        
        rej_h013not2[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
      }
      
      if( z23 >= func$upperbounds[2] & z2 >= func$upperbounds[1] &
          z1 < func$lowerbounds[1]){
        
        rej_h023not1[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
      }
      
      if( z23 < func$upperbounds[2] & z1 >= func$upperbounds[1] &
          z2 >=  func$upperbounds[1]){
        
        rej_h012not3[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
      }
      
      if( z23 < func$upperbounds[2] & z1 <= func$upperbounds[1] &
          z2 >=  func$upperbounds[1]){
        
        rej_h02not13[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
      }
      if( z23 < func$upperbounds[2] & z2 <= func$upperbounds[1] &
          z1 >=  func$upperbounds[1]){
        
        rej_h01not23[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
      }
      
      if( z23 < func$upperbounds[2] & z1 < func$lowerbounds[1] &
          z2 < func$lowerbounds[1]){
        
        bothineff_all[k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
      }
      
      
      
    }
    
    
  }
  
  
  # Average Probability of rejecting the null hypotheses
  
  prop_rej[[1]][k,1] <- sum(rej_h01not23[[1]][k,])/nsim
  
  prop_rej[[1]][k,2] <- sum(rej_h02not13[[1]][k,])/nsim
  
  prop_rej[[1]][k,3] <- sum(rej_h03not12[[1]][k,])/nsim
  
  prop_rej[[1]][k,4] <- sum(rej_h012not3[[1]][k,])/nsim
  
  prop_rej[[1]][k,5] <- sum(rej_h013not2[[1]][k,])/nsim
  
  prop_rej[[1]][k,6] <- sum(rej_h023not1[[1]][k,])/nsim
  
  prop_rej[[1]][k,7] <- sum(rej_h01andh02andh03[[1]][k,])/nsim
  
  prop_rej[[2]][k,1] <- sum(rej_h01not23[[2]][k,])/nsim
  
  prop_rej[[2]][k,2] <- sum(rej_h02not13[[2]][k,])/nsim
  
  prop_rej[[2]][k,3] <- sum(rej_h03not12[[2]][k,])/nsim
  
  prop_rej[[2]][k,4] <- sum(rej_h012not3[[2]][k,])/nsim
  
  prop_rej[[2]][k,5] <- sum(rej_h013not2[[2]][k,])/nsim
  
  prop_rej[[2]][k,6] <- sum(rej_h023not1[[2]][k,])/nsim
  
  prop_rej[[2]][k,7] <- sum(rej_h01andh02andh03[[2]][k,])/nsim
  
  prop_rej[[2]][k,8] <- sum(bothineff_all[k,])/nsim
  
  estimatedsamplesize[k] <- mean(estimatedss[k,])
  
  durationtrial[k] <- mean(dur[k,])
  
  timetofirstposclaim[k] <- mean(tfpc[k,])
}

colnames(prop_rej[[1]]) <- c(
  "Allprom_rejH01notH02H03_1stage", 
  "Allprom_rejH02notH01H03_1stage",
  "Allprom_rejH03notH01H02_1stage",
  "Allprom_rejH01andH02notH03_1stage", 
  "Allprom_rejH01andH03notH02_1stage",
  "Allprom_rejH02andH03notH01_1stage", 
  "Allprom_rejH01andH02andH03_1stage")

colnames(prop_rej[[2]]) <- c(
  "Allprom_rejH01notH02H03_2stage", 
  "Allprom_rejH02notH01H03_2stage",
  "Allprom_rejH03notH01H02_2stage",
  "Allprom_rejH01andH02notH03_2stage", 
  "Allprom_rejH01andH03notH02_2stage",
  "Allprom_rejH02andH03notH01_2stage", 
  "Allprom_rejH01andH02andH03_2stage",
  "Allprom_rejneitherH01andH02andH03")


colnames(estimatedsamplesize) <- "ESS Ordered ORD"

colnames(durationtrial) <- "Expected duration (months)"

colnames(timetofirstposclaim) <- "Expected TFPC (months)"

summary_MAMS <- cbind(popord,
                      sum(ceiling((rhopower*2*popord))),
                      apply(actualmaxss,1,mean, na.rm=TRUE),
                      scenario,
                      prop_rej, 
                      estimatedsamplesize,
                      durationtrial,
                      timetofirstposclaim,
                      func$upperbounds[1],
                      func$upperbounds[2],
                      func$lowerbounds[1])

colnames(summary_MAMS)[1:3]<- c("patientscontrol_perstage ORD",
                                "Max theoretical SS",
                                "Actual SS")

summary_MAMS <- summary_MAMS %>% mutate(
  
  
  
  Allprom_rejH01notH02H03 = Allprom_rejH01notH02H03_1stage+Allprom_rejH01notH02H03_2stage,
  
  Allprom_rejH02notH01H03 = Allprom_rejH02notH01H03_1stage+Allprom_rejH02notH01H03_2stage,
  
  Allprom_rejH03notH01H02 = Allprom_rejH03notH01H02_1stage+Allprom_rejH03notH01H02_2stage,
  
  Allprom_rejH01andH02notH03 = Allprom_rejH01andH02notH03_1stage+Allprom_rejH01andH02notH03_2stage,
  
  Allprom_rejH01andH03notH02 = Allprom_rejH01andH03notH02_1stage+Allprom_rejH01andH03notH02_2stage,
  
  Allprom_rejH02andH03notH01 = Allprom_rejH02andH03notH01_1stage+Allprom_rejH02andH03notH01_2stage,
  
  Allprom_rejH01andH02andH03 = Allprom_rejH01andH02andH03_1stage+Allprom_rejH01andH02andH03_2stage,
  
  Allprom_rejH01 = Allprom_rejH01notH02H03+Allprom_rejH01andH02notH03+Allprom_rejH01andH02andH03+Allprom_rejH01andH03notH02,
  
  Allprom_rejH02 = Allprom_rejH02notH01H03+Allprom_rejH02andH03notH01+Allprom_rejH01andH02andH03+Allprom_rejH01andH02notH03,
  
  Allprom_rejH03 = Allprom_rejH03notH01H02+Allprom_rejH02andH03notH01+Allprom_rejH01andH02andH03+Allprom_rejH01andH03notH02,
  
  Allprom_rejH01andH02 =  Allprom_rejH01andH02notH03+Allprom_rejH01andH02andH03,
  
  Allprom_rejH01andH03 =  Allprom_rejH01andH03notH02+Allprom_rejH01andH02andH03,
  
  Allprom_rejH02andH03 =  Allprom_rejH02andH03notH01+Allprom_rejH01andH02andH03,
  
  Allprom_FWER = Allprom_rejH01notH02H03+Allprom_rejH02notH01H03+Allprom_rejH03notH01H02+Allprom_rejH01andH02notH03+
    Allprom_rejH01andH03notH02+Allprom_rejH02andH03notH01+Allprom_rejH01andH02andH03,
  
  sumprob = Allprom_FWER+Allprom_rejneitherH01andH02andH03,
  
  `reject all` = Allprom_rejH01andH02andH03,
  
  `reject at least one` = Allprom_FWER,
  
  `reject L&M` = Allprom_rejH01andH02,
  
  u1 = func$upperbounds[1],
  u2 = func$upperbounds[2],
  l1 = func$lowerbounds[1]
  
)%>% mutate(
  
  p = paste("(", p0, ", ", p1,", ",p2,", ",p3,")",sep = ""),
  theta = paste("(", round(theta1, 3), ", ", round(theta2, 3),", ",round(theta3, 3),")",sep = "")
)%>% select(
  
  `patientscontrol_perstage ORD`,
  `Max theoretical SS`,
  `Actual SS`,
  `ESS Ordered ORD`,
  p,
  theta,
  `reject at least one`,
  `reject L&M`,
  `reject all`,
  Allprom_rejneitherH01andH02andH03,
  `Expected duration (months)`,
  `Expected TFPC (months)`,
  u1,
  u2,
  l1
  
  
)

write.csv(summary_MAMS, paste("4arm2stageMAMS_TFPC_new3_recrutrate", recrate,"_", Sys.Date(), 
                              ".csv", sep = ""),
          row.names = F)
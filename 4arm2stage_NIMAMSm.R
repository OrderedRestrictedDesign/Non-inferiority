# 4-arm 2-stage MAMS(m) design

rm(list=ls())

library(dplyr)


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

pop <- func$`sample size per arm per stage`#seq(func$`sample size per arm per stage`+7, func$`sample size per arm per stage`+8, by = 1)

upper <- func$a#seq(func$a+.03, func$a+.03, by = 0.001)

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
    
    func$upperbounds[1] <- ifelse(ushape == "obf", (upper[n] * 1/sqrt(ratio))[1],
                 ifelse(ushape == "pocock",rep(upper[n], stage-1),
                        (upper[n] * (1 + ratio/max(ratio))/sqrt(ratio))[1]))
    
    func$lowerbounds[1] <- ifelse(ushape == "obf", (-upper[n] * 1/sqrt(ratio)),
                 ifelse(ushape == "pocock",rep(-upper[n], stage-1),
                        (-upper[n] * ((1 -3* ratio/max(ratio))/sqrt(ratio)))[1]))
    
    func$upperbounds[2] <- ifelse(ushape == "obf", (upper[n] * 1/sqrt(ratio))[2],
                 ifelse(ushape == "pocock",rep(upper[n], stage-1),
                        (upper[n] * (1 + ratio/max(ratio))/sqrt(ratio))[2]))
    
    
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

estimatedss_allp <- matrix(data=NA, nrow=nrow(scenario),
                           ncol=nsim)

estimatedsamplesize_allp <- matrix(data=NA, ncol=1,
                                   nrow = nrow(scenario))

for (k in 1:nrow(scenario)){
  
  set.seed(seed = seed)
  
  mean0 <- scenario[k,1]
  
  mean1 <- scenario[k,2]
  
  mean2 <- scenario[k,3]
  
  mean3 <- scenario[k,4]
  
  for (j in 1:nsim){
    
    # Generate total population for each scenario and each simulation
    
    plac1st <-  rbinom(n = popord, size = 1, prob = mean0)
    
    dur1st <-  rbinom(n = popord, size = 1, prob = mean1)
    
    dur2st <-  rbinom(n = popord, size = 1, prob = mean2)
    
    dur3st <-  rbinom(n = popord, size = 1, prob = mean3)
    
    estimatedss_allp[k,j] <- length(plac1st)+length(dur1st)+length(dur2st)+length(dur3st)
    
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
    
    # Count when null hypothesis are rejected at the end of the stage
    
    
    if(
      (z1 >= func$upperbounds[1] & z2 >= func$upperbounds[1] &  z3 >= func$upperbounds[1])){ 
      
      
      rej_h01andh02andh03[[1]][k,j] = TRUE
      
    }
    
    
    if((z1 >= func$upperbounds[1] & z3 < func$lowerbounds[1] &  z2 < func$lowerbounds[1])){


      rej_h01not23[[1]][k,j] = TRUE

    }

    if((z2 >= func$upperbounds[1] & z1 < func$lowerbounds[1] &  z3 < func$lowerbounds[1])){
      
      
      rej_h02not13[[1]][k,j] = TRUE
      
    }
    
    if((z3 >= func$upperbounds[1] & z1 < func$lowerbounds[1] &  z2 < func$lowerbounds[1])){
      
      
      rej_h03not12[[1]][k,j] = TRUE
      
    }
    
    if((z1 >= func$upperbounds[1] & z3 < func$lowerbounds[1] &  z2 >= func$upperbounds[1])){
      
      
      rej_h012not3[[1]][k,j] = TRUE
      
    }
    
    if((z2 >= func$upperbounds[1] & z1 < func$lowerbounds[1] &  z3 >= func$upperbounds[1])){
      
      
      rej_h023not1[[1]][k,j] = TRUE
      
    }
    
    if((z1 >= func$upperbounds[1] & z2 < func$lowerbounds[1] &  z3 >= func$upperbounds[1])){
      
      
      rej_h013not2[[1]][k,j] = TRUE
      
    }
    
    if (z1 < func$lowerbounds[1] & z2 <func$lowerbounds[1] & z3 < func$lowerbounds[1]){
      
      
      bothineff_all[k,j] = TRUE
      
      
      
    }
    
    # Continue to the second stage with all arms
    
    if(z1 < func$upperbounds[1] & z1 > func$lowerbounds[1] &
       z2 < func$upperbounds[1] & z2 > func$lowerbounds[1] &
       z3 < func$upperbounds[1] & z3 > func$lowerbounds[1] ){
      
      # Estimated response rates
      
      plac21_2 <- rbinom(n = popord, size = 1, prob = mean0)
      
      dur21_2 <- rbinom(n = popord, size = 1, prob = mean1)
      
      dur22_2 <- rbinom(n = popord, size = 1, prob = mean2)
      
      dur23_2 <- rbinom(n = popord, size = 1, prob = mean3)
      
      estimatedss_allp[k,j] <-  estimatedss_allp[k,j]+length(plac21_2)+length(dur21_2)+length(dur22_2)+length(dur23_2)
      
      plac21 <- c(plac21_2, plac1st)
      
      dur21 <- c(dur21_2, dur1st)
      
      dur22 <- c(dur22_2, dur2st)
      
      dur23 <- c(dur23_2, dur3st)
      
      # Estimated response rates
      
      p_dur21 <- mean(dur21)
      
      p_plac21 <- mean(plac21)
      
      p_dur22 <- mean(dur22)
      
      p_dur23 <- mean(dur23)
      
      var21 <- (mean1*(1-mean1)/length(dur21))+(mean0*(1-mean0)/length(plac21))
      
      var22 <- (mean2*(1-mean2)/length(dur22))+(mean0*(1-mean0)/length(plac21))
      
      var23 <- (mean3*(1-mean3)/length(dur23))+(mean0*(1-mean0)/length(plac21))
      
      # Z-statistics for the two treatment durations
      
      z21 <- (p_dur21-p_plac21+noninf)/sqrt(var21)
      
      z22 <- (p_dur22-p_plac21+noninf)/sqrt(var22)
      
      z23 <- (p_dur23-p_plac21+noninf)/sqrt(var23)
      
      
      # Rejection of the null hypothesis at the second stage
      
      
      
      
      if((z22 < func$upperbounds[2] & z23 < func$upperbounds[2] &  z21 >= func$upperbounds[2])) {
        
        rej_h01not23[[2]][k,j] = TRUE
      }
      
      if(
        (z22 >= func$upperbounds[2] & z23 < func$upperbounds[2] &  z21 >= func$upperbounds[2])){
        
        rej_h012not3[[2]][k,j] = TRUE
        
      } 
      if(
        (z22 < func$upperbounds[2] & z23 >= func$upperbounds[2] &  z21 >= func$upperbounds[2])){
        
        rej_h013not2[[2]][k,j] = TRUE
      }
      
      if(
        (z21 >= func$upperbounds[2] & z22 >= func$upperbounds[2] &  z23 >= func$upperbounds[2])
      ){ 
        
        
        rej_h01andh02andh03[[2]][k,j] = TRUE
        
        
      }
      
      
      if((z21 < func$upperbounds[2] & z23 < func$upperbounds[2] &  z22 >= func$upperbounds[2])){
        
        rej_h02not13[[2]][k,j] = TRUE
      }
      
      if( 
        (z21 < func$upperbounds[2] & z23 >= func$upperbounds[2] &  z22 >= func$upperbounds[2]) ){
        
        rej_h023not1[[2]][k,j] = TRUE
      }
      
      
      # rej H03
      
      
      if((z21 < func$upperbounds[2] & z22 < func$upperbounds[2] &  z23 >= func$upperbounds[2])){
        
        rej_h03not12[[2]][k,j] = TRUE
      }
      
      
      
      
      
      if (z21 < func$upperbounds[2] & z22 <func$upperbounds[2] & z23 < func$upperbounds[2]){
        
        
        bothineff_all[k,j] = TRUE
        
        
        
      }
      
      
      
    }
    
    
    # # Continue to the second stage with first and second arm
    
    if((z1 < func$upperbounds[1] & z1 >=  func$lowerbounds[1] &
       z2 < func$upperbounds[1] & z2 >=  func$lowerbounds[1] &
       z3 >=  func$upperbounds[1]) ||
       (z1 < func$upperbounds[1] & z1 >=  func$lowerbounds[1] &
        z2 < func$upperbounds[1] & z2 >=  func$lowerbounds[1] &
        z3 <=  func$lowerbounds[1])){
      
      # Estimated response rates
      
      plac21_2 <- rbinom(n = popord, size = 1, prob = mean0)
      
      dur21_2 <- rbinom(n = popord, size = 1, prob = mean1)
      
      dur22_2 <- rbinom(n = popord, size = 1, prob = mean2)
      
      estimatedss_allp[k,j] <-  estimatedss_allp[k,j]+length(plac21_2)+length(dur21_2)+length(dur22_2)
      
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
      
      if( ( z21 >= func$upperbounds[2] & z22 < func$upperbounds[2] & z3 >= func$upperbounds[1]) ){
        
        rej_h013not2[[2]][k,j] = TRUE
        
        
        
      }
      
      if( ( z21 >= func$upperbounds[2] & z22 < func$upperbounds[2] & z3 < func$lowerbounds[1]) ){
        
        rej_h01not23[[2]][k,j] = TRUE
        
        
        
      }
      
      if( ( z22 >= func$upperbounds[2] & z21 < func$upperbounds[2] & z3 >= func$upperbounds[1]) ){
        
        rej_h023not1[[2]][k,j] = TRUE
        
        
        
      }
      
      if( ( z22 >= func$upperbounds[2] & z21 < func$upperbounds[2] & z3< func$lowerbounds[1]) ){
        
        rej_h02not13[[2]][k,j] = TRUE
        
        
        
      }
      
      if(( z21 >= func$upperbounds[2] & z22 >= func$upperbounds[2]) & z3 >= func$upperbounds[1]){
        
        rej_h01andh02andh03[[2]][k,j] = TRUE
        
        
      }
      
      if(( z21 >= func$upperbounds[2] & z22 >= func$upperbounds[2]) & z3 < func$lowerbounds[1]){
        
        rej_h012not3[[2]][k,j] = TRUE
        
        
      }
      
      
      if( ( z21 < func$upperbounds[2] & z22 < func$upperbounds[2]) & z3 >= func$upperbounds[1]){
        
        rej_h03not12[[2]][k,j] = TRUE
        
        
        
      }
      
      if( ( z21 < func$upperbounds[2] & z22 < func$upperbounds[2]) & z3 < func$lowerbounds[1]){
        
        bothineff_all[k,j] = TRUE
        
        
        
      }
      
      
      
      
      
    }
    
    
    # Continue to the second stage with first and third arm
    
    if((z1 < func$upperbounds[1] & z1 >=  func$lowerbounds[1] &
       z3 < func$upperbounds[1] & z3 >=  func$lowerbounds[1] &
       z2 >=  func$upperbounds[1]) ||
       (z1 < func$upperbounds[1] & z1 >=  func$lowerbounds[1] &
        z3 < func$upperbounds[1] & z3 >=  func$lowerbounds[1] &
        z2 <  func$lowerbounds[1])){
      
      # Estimated response rates
      
      plac21_2 <- rbinom(n = popord, size = 1, prob = mean0)
      
      dur21_2 <- rbinom(n = popord, size = 1, prob = mean1)
      
      dur23_2 <- rbinom(n = popord, size = 1, prob = mean3)
      
      estimatedss_allp[k,j] <-  estimatedss_allp[k,j]+length(plac21_2)+length(dur21_2)+length(dur23_2)
      
      plac21 <- c(plac21_2, plac1st)
      
      dur21 <- c(dur21_2, dur1st)

      dur23 <- c(dur23_2, dur3st)
      
      # Estimated response rates
      
      p_dur21 <- mean(dur21)
      
      p_plac21 <- mean(plac21)
      
      p_dur23 <- mean(dur23)
      
      var21 <- (mean1*(1-mean1)/length(dur21))+(mean0*(1-mean0)/length(plac21))

      var23 <- (mean3*(1-mean3)/length(dur23))+(mean0*(1-mean0)/length(plac21))
      
      # Z-statistics for the two treatment durations
      
      z21 <- (p_dur21-p_plac21+noninf)/sqrt(var21)
      
      z23 <- (p_dur23-p_plac21+noninf)/sqrt(var23)
      
      
      if( ( z21 >= func$upperbounds[2] & z23 < func$upperbounds[2] & z2 >=  func$upperbounds[1]) ){
        
        rej_h012not3[[2]][k,j] = TRUE
        
        
        
      }
      
      if( ( z21 >= func$upperbounds[2] & z23 < func$upperbounds[2] & z2 <  func$lowerbounds[1]) ){
        
        rej_h01not23[[2]][k,j] = TRUE
        
        
        
      }
      
      if( ( z23 >= func$upperbounds[2] & z21 < func$upperbounds[2] & z2 >=  func$upperbounds[1]) ){
        
        rej_h023not1[[2]][k,j] = TRUE
        
        
        
      }
      
      if( ( z23 >= func$upperbounds[2] & z21 < func$upperbounds[2] & z2 <  func$lowerbounds[1]) ){
        
        rej_h03not12[[2]][k,j] = TRUE
        
        
        
      }
      
      if(
        ( z21 >= func$upperbounds[2] & z23 >= func$upperbounds[2]  & z2 >=  func$upperbounds[1])){
        
        
        rej_h01andh02andh03[[2]][k,j] = TRUE
        
        
      }
      
      if(
        ( z21 >= func$upperbounds[2] & z23 >= func$upperbounds[2]  & z2<  func$lowerbounds[1])){
        
        
        rej_h013not2[[2]][k,j] = TRUE
        
        
      }
      
      
      if( ( z21 < func$upperbounds[2] & z23 < func$upperbounds[2]& z2 >=  func$upperbounds[1]) ){
        
        rej_h02not13[[2]][k,j] = TRUE
        
        
        
      }
      
      if( ( z21 < func$upperbounds[2] & z23 < func$upperbounds[2]& z2 <  func$lowerbounds[1]) ){
        
        bothineff_all[k,j] = TRUE
        
        
        
      }
      
      
    }
    
    # Continue to the second stage with third and second arm
    
    if((z3 < func$upperbounds[1] & z3 >=  func$lowerbounds[1] &
       z2 < func$upperbounds[1] & z2 >=  func$lowerbounds[1] &
       z1 >=  func$upperbounds[1])||
       (z3 < func$upperbounds[1] & z3 >=  func$lowerbounds[1] &
        z2 < func$upperbounds[1] & z2 >=  func$lowerbounds[1] &
        z1 <  func$lowerbounds[1])){
      
      # Estimated response rates
      
      plac21_2 <- rbinom(n = popord, size = 1, prob = mean0)

      dur22_2 <- rbinom(n = popord, size = 1, prob = mean2)
      
      dur23_2 <- rbinom(n = popord, size = 1, prob = mean3)
      
      estimatedss_allp[k,j] <-  estimatedss_allp[k,j]+length(plac21_2)+length(dur22_2)+length(dur23_2)
      
      plac21 <- c(plac21_2, plac1st)
      
      dur22 <- c(dur22_2, dur2st)
      
      dur23 <- c(dur23_2, dur3st)
      
      # Estimated response rates

      p_plac21 <- mean(plac21)
      
      p_dur22 <- mean(dur22)
      
      p_dur23 <- mean(dur23)

      var22 <- (mean2*(1-mean2)/length(dur22))+(mean0*(1-mean0)/length(plac21))
      
      var23 <- (mean3*(1-mean3)/length(dur23))+(mean0*(1-mean0)/length(plac21))
      
      # Z-statistics for the two treatment durations
      

      z22 <- (p_dur22-p_plac21+noninf)/sqrt(var22)
      
      z23 <- (p_dur23-p_plac21+noninf)/sqrt(var23)
      
      
      if( ( z22 >= func$upperbounds[2] & z23 < func$upperbounds[2] &
            z1 >=  func$upperbounds[1]) ){
        
        rej_h012not3[[2]][k,j] = TRUE
        
        
        
      }
      
      if( ( z22 >= func$upperbounds[2] & z23 < func$upperbounds[2] &
            z1 <  func$lowerbounds[1]) ){
        
        rej_h023not1[[2]][k,j] = TRUE
        
        
        
      }
      
      if( ( z23 >=func$upperbounds[2] & z22 < func$upperbounds[2] &
            z1 >=  func$upperbounds[1])  ){
        
        rej_h013not2[[2]][k,j] = TRUE
        
        
        
      }
      
      if( ( z23 >=func$upperbounds[2] & z22 < func$upperbounds[2] &
            z1 <  func$lowerbounds[1])){
        
        rej_h03not12[[2]][k,j] = TRUE
        
        
        
      }
      
      if( 
        ( z22 >= func$upperbounds[2] & z23 >= func$upperbounds[2]&
          z1 >=  func$upperbounds[1])){
        
        
        rej_h01andh02andh03[[2]][k,j] = TRUE
        
        
      }
      
      if( 
        ( z22 >= func$upperbounds[2] & z23 >= func$upperbounds[2]&
          z1 <  func$lowerbounds[1])){
        
        
        rej_h023not1[[2]][k,j] = TRUE
        
        
      }
      
      
      if( ( z22 < func$upperbounds[2] & z23 < func$upperbounds[2]&
            z1 >=  func$upperbounds[1]) ){
        
        rej_h01not23[[2]][k,j] = TRUE
        
        
        
      }
      
      if( ( z22 < func$upperbounds[2] & z23 < func$upperbounds[2]&
            z1<  func$lowerbounds[1]) ){
        
        bothineff_all[k,j] = TRUE
        
        
        
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
      
      # Estimated response rates
      
      plac21_2 <- rbinom(n = popord, size = 1, prob = mean0)
      
      dur21_2 <- rbinom(n = popord, size = 1, prob = mean1)
     
      estimatedss_allp[k,j] <-  estimatedss_allp[k,j]+length(plac21_2)+length(dur21_2)
      
      plac21 <- c(plac21_2, plac1st)
      
      dur21 <- c(dur21_2, dur1st)
      
       # Estimated response rates
      
      p_dur21 <- mean(dur21)
      
      p_plac21 <- mean(plac21)
      
      
      var21 <- (mean1*(1-mean1)/length(dur21))+(mean0*(1-mean0)/length(plac21))
      
      
      # Z-statistics for the two treatment durations
      
      z21 <- (p_dur21-p_plac21+noninf)/sqrt(var21)
      
     
      # rej H01
      
      if( z21 >= func$upperbounds[2] & z2 >= func$upperbounds[1] &
          z3 >=  func$upperbounds[1] ){
        
        
        rej_h01andh02andh03[[2]][k,j] = TRUE
      }
      
      if( z21 >= func$upperbounds[2] & z2 < func$lowerbounds[1] &
          z3 <  func$lowerbounds[1] ){
        
        
        rej_h01not23[[2]][k,j] = TRUE
      }
      
      if( z21 >= func$upperbounds[2] & z2 >= func$upperbounds[1] &
          z3 <  func$lowerbounds[1] ){
        
        
        rej_h012not3[[2]][k,j] = TRUE
      }
      
      if( z21 >= func$upperbounds[2] & z3 >= func$upperbounds[1] &
          z2 <  func$lowerbounds[1] ){
        
        
        rej_h013not2[[2]][k,j] = TRUE
      }
      
      if( z21 < func$upperbounds[2] & z2 >= func$upperbounds[1] &
          z3 >=  func$upperbounds[1]){
        
        
        rej_h023not1[[2]][k,j] = TRUE
      }
      
      if( z21 < func$upperbounds[2] & z2 >= func$upperbounds[1] &
          z3 <  func$lowerbounds[1]){
        
        
        rej_h02not13[[2]][k,j] = TRUE
      }
      
      if( z21 < func$upperbounds[2] & z3 >= func$upperbounds[1] &
          z2 <  func$lowerbounds[1]){
        
        
        rej_h03not12[[2]][k,j] = TRUE
      }
      
      if( z21 < func$upperbounds[2] & z2 < func$lowerbounds[1] &
          z3 <  func$lowerbounds[1] ){
        
        
        bothineff_all[k,j] = TRUE
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
      
      # Estimated response rates
      
      plac21_2 <- rbinom(n = popord, size = 1, prob = mean0)

      dur22_2 <- rbinom(n = popord, size = 1, prob = mean2)
      
      estimatedss_allp[k,j] <-  estimatedss_allp[k,j]+length(plac21_2)+length(dur22_2)
      
      plac21 <- c(plac21_2, plac1st)
      
      dur22 <- c(dur22_2, dur2st)
      
      # Estimated response rates
      
      
      p_plac21 <- mean(plac21)
      
      p_dur22 <- mean(dur22)
      
      
      var22 <- (mean2*(1-mean2)/length(dur22))+(mean0*(1-mean0)/length(plac21))
      
      # Z-statistics for the two treatment durations
      
      
      z22 <- (p_dur22-p_plac21+noninf)/sqrt(var22)
      
      # rej H02
      
      if( z22 >= func$upperbounds[2] & z1 >= func$upperbounds[1] &
          z3 >=  func$upperbounds[1]){
        
        rej_h01andh02andh03[[2]][k,j] = TRUE
        
      }
      
      if( z22 >= func$upperbounds[2] & z1 < func$lowerbounds[1] &
          z3 < func$lowerbounds[1]){
        
        rej_h02not13[[2]][k,j] = TRUE
        
      }
      
      if( z22 >= func$upperbounds[2] & z1 >= func$upperbounds[1] &
          z3 < func$lowerbounds[1]){
        
        rej_h012not3[[2]][k,j] = TRUE
        
      }
      
      if( z22 >= func$upperbounds[2] & z3 >= func$upperbounds[1] &
          z1 < func$lowerbounds[1]){
        
        rej_h023not1[[2]][k,j] = TRUE
        
      }
      
      if( z22 < func$upperbounds[2]& z1 >= func$upperbounds[1] &
          z3 >=  func$upperbounds[1]){
        
        rej_h013not2[[2]][k,j] = TRUE
        
      }
      
      if( z22 < func$upperbounds[2]& z1 < func$lowerbounds[1] &
          z3 >=  func$upperbounds[1]){
        
        rej_h03not12[[2]][k,j] = TRUE
        
      }
      
      if( z22 < func$upperbounds[2]& z3 < func$lowerbounds[1] &
          z1 >=  func$upperbounds[1]){
        
        rej_h01not23[[2]][k,j] = TRUE
        
      }
      
      if( z22 < func$upperbounds[2]& z1 < func$lowerbounds[1] &
          z3 < func$lowerbounds[1]){
        
        bothineff_all[k,j] = TRUE
        
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
      
      # Estimated response rates
      
      plac21_2 <- rbinom(n = popord, size = 1, prob = mean0)
      
      dur23_2 <- rbinom(n = popord, size = 1, prob = mean3)
      
      estimatedss_allp[k,j] <-  estimatedss_allp[k,j]+length(plac21_2)+length(dur23_2)
      
      plac21 <- c(plac21_2, plac1st)
      
      dur23 <- c(dur23_2, dur3st)
      
      # Estimated response rates
      
      p_plac21 <- mean(plac21)
      
      p_dur23 <- mean(dur23)
      
      var23 <- (mean3*(1-mean3)/length(dur23))+(mean0*(1-mean0)/length(plac21))
      
      # Z-statistics for the two treatment durations

      z23 <- (p_dur23-p_plac21+noninf)/sqrt(var23)
      
      
      # rej H03
      
      if( z23 >= func$upperbounds[2] & z1 >= func$upperbounds[1] &
          z2 >=  func$upperbounds[1]){
        
        rej_h01andh02andh03[[2]][k,j] = TRUE
        
      }
      
      if( z23 >= func$upperbounds[2] & z1 < func$lowerbounds[1] &
          z2 < func$lowerbounds[1]){
        
        rej_h03not12[[2]][k,j] = TRUE
        
      }
      
      if( z23 >= func$upperbounds[2] & z1 >= func$upperbounds[1] &
          z2 < func$lowerbounds[1]){
        
        rej_h013not2[[2]][k,j] = TRUE
        
      }
      
      if( z23 >= func$upperbounds[2] & z2 >= func$upperbounds[1] &
          z1 < func$lowerbounds[1]){
        
        rej_h023not1[[2]][k,j] = TRUE
        
      }
      
      if( z23 < func$upperbounds[2] & z1 >= func$upperbounds[1] &
          z2 >=  func$upperbounds[1]){
        
        rej_h012not3[[2]][k,j] = TRUE
        
      }
      
      if( z23 < func$upperbounds[2] & z1 <= func$upperbounds[1] &
          z2 >=  func$upperbounds[1]){
        
        rej_h02not13[[2]][k,j] = TRUE
        
      }
      if( z23 < func$upperbounds[2] & z2 <= func$upperbounds[1] &
          z1 >=  func$upperbounds[1]){
        
        rej_h01not23[[2]][k,j] = TRUE
        
      }
      
      if( z23 < func$upperbounds[2] & z1 < func$lowerbounds[1] &
          z2 < func$lowerbounds[1]){
        
        bothineff_all[k,j] = TRUE
        
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
  
  estimatedsamplesize_allp[k] <- mean(estimatedss_allp[k,]) 
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


colnames(estimatedsamplesize_allp) <- "ESS All promising"

summary_MAMS <- cbind(popord,popord*2*arms,scenario,prop_rej, estimatedsamplesize_allp,
                     func$upperbounds[1],
                     func$upperbounds[2],
                     func$lowerbounds[1])

colnames(summary_MAMS)[1:2] <- c("patients paps", "total sample")

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
  
) %>% select(
  
  `patients paps`,
  
  `total sample`,
  
  `ESS All promising`,
  
  theta1,
  
  theta2,
  
  theta3,
  
  Allprom_rejH01,
  
  Allprom_rejH02,
  
  Allprom_rejH03,
  
  Allprom_rejH01andH02,
  
  Allprom_rejH01andH03,
  
  Allprom_rejH02andH03,
  
  Allprom_rejH01andH02andH03,
  
  Allprom_rejH01notH02H03,
  
  Allprom_rejH02notH01H03,
  
  Allprom_rejH03notH01H02,
  
  Allprom_rejH01andH02notH03,
  
  Allprom_rejH01andH03notH02,
  
  Allprom_rejH02andH03notH01,

  Allprom_FWER,
  
  Allprom_rejneitherH01andH02andH03,
  
  sumprob,
  
  `reject all`,
  
  `reject at least one`,
  
  `reject L&M`,
  
  u1,
  
  u2,
  
  l1
  
)

typeI <- summary_MAMS[1,]$`Allprom_FWER`

power <- summary_MAMS[2,powereq]

if(typeI <= alpha & power > 1-beta ){
  
  mat2 <- rbind(mat2, cbind(pop[u],
                            upper[n],
                            func$upperbounds[1],
                            func$upperbounds[2],
                            func$lowerbounds[1],
                            typeI,
                            power))
  
 
  
}



  }
  
  
}

colnames(mat2) <- c("sample size", "a","u1", "u2", "l1", "typeI under null", "power")

mat2 <- cbind(mat2, pc,noninf)


rm(list=ls())

##### Time to first positive claim

####### Call for functions ####

##### Function to find critical bounds for 3-arm 2-stage MAMS(m) NI design ##############

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
#' @lfix: fixed value for lower bound

bounds_3arm2stagem <- function(u1, interimN, arms,r,
                               alpha, cov, first, ushape,lshape,lfix, prec){
  
  library(mvtnorm)
  library(gtools)
  
  alphat <- rep(0,length = interimN)
  
  alphatcaseup<- NULL
  
  alphatcasedown<- NULL
  
  alphatcasefin<- NULL
  
  for(j in 1:length(u1)){
    
    for (i in 1:interimN){
      
      
      L <- array(0, dim = i*(arms-1))
      U <- array(0, dim = i*(arms-1))
      
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
                               sigma = first,
                               algorithm = GenzBretz(abseps = 1*10^-12))[1]
        
      }
      if(i>1){
        
        
        alphat[i] <- pmvnorm(lower = c(l[1],-Inf,  u[2], -Inf),
                             upper = c(u[1],u[1],  Inf, Inf),
                             sigma = cov,
                             algorithm = GenzBretz(abseps = 1*10^-12))[1]+
          pmvnorm(lower = c(-Inf,l[1],  -Inf,u[2]),
                  upper = c(u[1],u[1], Inf, Inf),
                  sigma = cov,
                  algorithm = GenzBretz(abseps = 1*10^-12))[1]-
          pmvnorm(lower = c(l[1], l[1],u[2], u[2]),
                  upper = c(u[1],u[1], Inf, Inf),
                  sigma = cov,
                  algorithm = GenzBretz(abseps = 1*10^-12))[1]
        
        
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

##### Function to find sample size for 3-arm 2-stage MAMS(m) NI design ###############

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

boundaries_3arm2stagem <- function(theta,
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
                                   power){ 
  
  library(mvtnorm)
  library(cubature)
  
  interimN <- stage
  
  #construction of the covariance matrix
  
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
  
  
  #covariance matrix at the first stage
  
  firststagematrix <- covbfin[1:(arms-1),1:(arms-1)]
  
  diag(firststagematrix) <- 1
  
  ## Research of the critical bounds
  
  u1 <- seq(from = 7, to = 0, by = -1)
  
  first <- bounds_3arm2stagem(u1, stage, arms,r, 
                              alpha, covbfin, firststagematrix,
                              ushape, lshape,lfix,prec)
  
  low <- first$alow
  
  up <- first$aup
  
  p <- rep(1, times = prec)
  
  for (p in 1:length(p)){
    
    callf <- bounds_3arm2stagem(seq(from = low, 
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
  
  # Search of the sample size to reach the desired power
  
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
  
  
  pop <- seq(from = 1, to = 1500, by = 1)
  
  
  for (n in 1:length(pop)){
    
    sigma1 <- (varpc/(rhopower[3]*pop[n]))+(vart1/(rhopower[1]*pop[n]))
    
    sigma2 <- (varpc/(rhopower[3]*pop[n]))+(vart2/(rhopower[2]*pop[n]))
    
    sigman <- c(sigma1,sigma2)*rep(1/r, each=arms-1)
    
    mean2<- rep(((theta[2:(arms)]-theta[1])+noninf),2)/(sqrt(sigman))
    
    mean1 <- mean2[1:(arms-1)]
    
    for (i in 1:interimN){
      
      L <- array(0, dim = i*(arms-1))
      U <- array(0, dim = i*(arms-1))
      
      if (ushape == "obf") {
        u <- a * 1/sqrt(r)
      }
      else if (ushape == "pocock") {
        u <- rep(a, i)
      }
      else if (ushape == "triangular") {
        u <- a * (1 + r/max(r))/sqrt(r)
      }
      
      if (lshape == "obf") {
        l <- c(-a * 1/sqrt(r))[1:(i-1)]
      }
      else if (lshape == "pocock") {
        l <- c(rep(-a, i-1 ))
      }
      else if (lshape == "triangular") {
        if (ushape == "triangular") {
          l <- c(-a * ((1 -3* r/max(r))/sqrt(r)))[1:(i-1)]
        }
      }
      else if (lshape == "fixed") {
        l <- c(rep(lfix, i-1 ))
      }
      
      
      # For 3-arm 2-stage 
      
      if(power== "reject all"){
        
        if(i ==1){
          set.seed(123456)
          
          betat[i] <- pmvnorm(lower = rep(upperbound[1],times = arms-1),
                              upper = rep(Inf,times = arms-1),
                              sigma = firststagematrixB,
                              mean = mean1,
                              algorithm = GenzBretz(abseps = 1*10^-12))[1]
          
          
        }
        else{
          
          betat[i] <- pmvnorm(lower = c(l[1],u[1],  u[2], -Inf),
                              upper = c(u[1],Inf,  Inf, Inf),
                              sigma = covbfinB,
                              mean = mean2,
                              algorithm = GenzBretz(abseps = 1*10^-12))[1]+
            pmvnorm(lower = c(u[1],l[1],  -Inf,u[2]),
                    upper = c(Inf,u[1], Inf, Inf),
                    sigma = covbfinB,
                    mean = mean2,
                    algorithm = GenzBretz(abseps = 1*10^-12))[1]+
            pmvnorm(lower = c(l[1], l[1],u[2], u[2]),
                    upper = c(u[1],u[1], Inf, Inf),
                    sigma = covbfinB,
                    mean = mean2,
                    algorithm = GenzBretz(abseps = 1*10^-12))[1]
          
        }
        
        
      }
      if(power=="reject at least one"){
        
        if(i ==1){
          set.seed(123456)
          
          betat[i] <- 1-pmvnorm(lower = rep(-Inf,times = arms-1),
                                upper = rep(u[1],times = arms-1),
                                sigma = firststagematrixB,
                                mean = mean1,
                                algorithm = GenzBretz(abseps = 1*10^-12))[1]
          
        }
        if(i>1){
          
          
          
          betat[i] <- pmvnorm(lower = c(l[1],-Inf,  u[2], -Inf),
                              upper = c(u[1],u[1],  Inf, Inf),
                              sigma = covbfinB,
                              mean = mean2,
                              algorithm = GenzBretz(abseps = 1*10^-12))[1]+
            pmvnorm(lower = c(-Inf,l[1],  -Inf,u[2]),
                    upper = c(u[1],u[1], Inf, Inf),
                    sigma = covbfinB,
                    mean = mean2,
                    algorithm = GenzBretz(abseps = 1*10^-12))[1]-
            pmvnorm(lower = c(l[1], l[1],u[2], u[2]),
                    upper = c(u[1],u[1], Inf, Inf),
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
                  lowerbound,upperbound,maxsample, a )
  
  names(results) <- c("sum_alphat", 
                      "lowerbounds", "upperbounds", "sample size per arm per stage", "a")
  
  return(results)
}




#### SIMULATION 3 arms - 2 stages MAMS(m) NI design ######

# PARAMETERS for both designs

seed <- 64736 # set seed to simulate the data

nsim <- 5*10^4 # number of simulations

arms <- 3 # number of arms

stage <- 2 # number of stages

alpha <- 0.05 # alpha-level of the test

noninf <- 0.1 # non-inferiority margin

pc <- 0.86 # response rate on the cntrol arm

pt1 <-pc # response rate on the first treatment arm

pt2 <-pc # response rate on the second treatment arm

beta <- 0.2# beta-level

ushape <- "triangular" # shape of upper critical boundary

lshape <- "triangular" # shape of lower critical boundary

thetapow <- 0 # value of treatment effect for power configuration

powerconf <- c(0, rep(thetapow, 2)) # power configuration

rho <- c(1,1,1) #sample size in each treatment divided by sample size in the control: rho[1]=n11/n01, rho[2]=n21/n01, rho[3]=1

r <- ratio <- c(1,2) #allocation ratio sample size at different stages: n1=r[1], n2=r[2]*n1,...

prec <- 3 # precision - number of digits after decimal point - for the critical boundaries

powereq <- "reject at least one" # type of power 

# call function to find the sample size per arm per stage and the critical boundaries

func <- boundaries_3arm2stagem(theta = c(pc,pt1,pt2)+powerconf,
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
                               power = powereq)

func

# simulation scenario 

####### RESTRUCTURE scenarios

plac <- rep(pc, 9)

trt1 <- c(rep(pc, 8), pc-noninf)

trt2 <-  c(seq(pc-noninf, pc+noninf, by = 0.03),pc,pc-noninf)

scen <- cbind(plac,trt1,trt2)

scenario <- matrix(data = scen,
                   nrow = nrow(scen),
                   ncol = ncol(scen))

scenario <- cbind(scenario, scenario[,2]-scenario[,1], scenario[,3]-scenario[,1])

colnames(scenario) <- c("p0", "pL", "pS","thetaL", "thetaS")

scenario <- as.data.frame(scenario)

# probability to allocate every patient on each treatment arm

p <- rep(1/arms, arms)

recrate <- 30

# Durations of the treatment

durations <- c(6,4,3)

######## Sample size and bound grid values 

popord <- 74#func$`sample size per arm per stage`+1#seq(func$`sample size per arm per stage`-5, func$`sample size per arm per stage`+5, by = 1)

upper <- func$a#seq(func$a-.05, func$a+.05, by = 0.001)

## simulations

# critical boundaries

func$upperbounds[1] <- u1 <- ifelse(ushape == "obf", (upper * 1/sqrt(ratio))[1],
                                    ifelse(ushape == "pocock",rep(upper, stage-1),
                                           (upper * (1 + ratio/max(ratio))/sqrt(ratio))[1]))

func$lowerbounds[1] <- l1 <- ifelse(ushape == "obf", (-upper * 1/sqrt(ratio)),
                                    ifelse(ushape == "pocock",rep(-upper, stage-1),
                                           (-upper * ((1 -3* ratio/max(ratio))/sqrt(ratio)))[1]))

func$upperbounds[2] <- u2 <- ifelse(ushape == "obf", (upper * 1/sqrt(ratio))[2],
                                    ifelse(ushape == "pocock",rep(upper, stage-1),
                                           (upper * (1 + ratio/max(ratio))/sqrt(ratio))[2]))

summary_MAMS <- NULL



# Save rejections of H01 and H02

rej_h01andh02 <- list()

rej_h01andh02[[1]] <- matrix(data=FALSE, nrow=nrow(scenario),
                             ncol=nsim)

rej_h01andh02[[2]] <- matrix(data=FALSE, nrow=nrow(scenario),
                             ncol=nsim)

# Save rejections of H01 and not H02

rej_h01noh02 <- list()

rej_h01noh02[[1]] <- matrix(data=FALSE, nrow=nrow(scenario),
                            ncol=nsim)

rej_h01noh02[[2]] <- matrix(data=FALSE, nrow=nrow(scenario),
                            ncol=nsim)

# Save rejections of H02 and not H01

rej_h02noh01 <- list()

rej_h02noh01[[1]] <- matrix(data=FALSE, nrow=nrow(scenario),
                            ncol=nsim)

rej_h02noh01[[2]] <- matrix(data=FALSE, nrow=nrow(scenario),
                            ncol=nsim)

bothineff_all  <- matrix(data=FALSE, nrow=nrow(scenario),
                         ncol=nsim)

prop_rej <- list()

prop_rej[[1]] <- matrix(data=NA, ncol=3,
                        nrow = nrow(scenario))
prop_rej[[2]] <- matrix(data=NA, ncol=4,
                        nrow = nrow(scenario))

# Estimated sample size for each scenario. --> The average of the sample size used for each simulation

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

timetofirstposclaim <- matrix(data=NA, ncol=1,
                              nrow = nrow(scenario))

tfpc <- matrix(data=0, nrow=nrow(scenario),
               ncol=nsim)


for (k in 1:nrow(scenario)){
  
  set.seed(seed = seed)
  
  mean0 <- scenario[k,1]
  
  mean1 <- scenario[k,2]
  
  mean2 <- scenario[k,3]
  
  for (j in 1:nsim){
    
    plac1st <- dur1st <- dur2st <- NULL
    
    plac21 <- dur21 <- dur22 <- NULL
    
    plac_2missing <- short <- long <- plac <- 0
    
    count <- NULL
    
    lastpatplac <- FALSE # when I recruit the last patient on the placebo
    
    lastpatlong <- FALSE # when I recruit the last patient on the longest
    
    lastpatshort <- FALSE # when I recruit the last patient on the shortest
    
    #recruit until get popord on the shortest
    
    while(plac < popord || short < popord || long < popord){
      
      count <- c(count,sample(c("C","L","S"), recrate, prob = p, replace=TRUE))
      short <- sum(count=="S")
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
      
      
    }
    #table(count)
    
    n_plac <- sum(count=="C")
    n_long <- sum(count=="L")
    n_short <- sum(count=="S")
    months <- length(count)/recrate
    
    n_interim <- min(c(n_plac,n_long,n_short))
    
    plac11st <- rbinom(n = n_interim, prob= mean0, size = 1)
    
    dur11st <- rbinom(n = n_interim, prob= mean1, size = 1)
    
    dur21st <- rbinom(n = n_interim, prob= mean2, size = 1)
    
    # continue to recruit until observe the patients on the shortest duration
    
    plac_over <- n_plac
    
    count_over <- NULL
    
    lastpatientsobserved <- c(month_lastpatplac, month_lastpatlong,month_lastpatshort)
    
    time_interim_month <- max(c(month_lastpatplac+durations[1],#+6,
                                month_lastpatlong+durations[2],#+6,
                                month_lastpatshort+durations[3]))#+6))
    
    index_timeinter <- which.max(c(month_lastpatplac+durations[1],#+6,
                                   month_lastpatlong+durations[2],#+6,
                                   month_lastpatshort+durations[3]))#+6))
    
    continuerecr <- time_interim_month-lastpatientsobserved[index_timeinter]
    
    dur[k,j] <- time_interim_month
    
    #continue recruitment while waiting for interim
    
    while((length(count_over)/recrate) < continuerecr & plac_over < 2*popord){
      count_over <- c(count_over,sample(c("C","L","S"), recrate, prob = p, replace=TRUE))
      months_over <- length(count_over)/recrate
      plac_over <- sum(count_over=="C")
    }
    #table(count_over)
    
    n_plac_over <- sum(count_over=="C")
    n_long_over <- sum(count_over=="L")
    n_short_over <- sum(count_over=="S")
    
    plac1st <- c(plac11st, rbinom(n = (n_plac-n_interim)+n_plac_over, prob= mean0, size = 1))
    
    dur1st <- c(dur11st, rbinom(n = (n_long-n_interim)+n_long_over, prob= mean1, size = 1))
    
    dur2st <- c(dur21st, rbinom(n = (n_short-n_interim)+n_short_over, prob= mean2, size = 1))
    
    estimatedss[k,j] <- length(plac1st)+length(dur1st)+length(dur2st)
    
    #dur[k,j] <- dur[k,j] +  days_over
    
    # Estimated response rates on the observed data
    
    p_plac11 <- mean(plac11st)
    
    p_dur11 <- mean(dur11st)
    
    p_dur12 <- mean(dur21st)
    
    
    p1 <- (mean1*(1-mean1)/length(dur11st))+(mean0*(1-mean0)/length(plac11st))
    
    p2 <- (mean2*(1-mean2)/length(dur21st))+(mean0*(1-mean0)/length(plac11st))
    
    sdpooled11 <- p1
    
    sdpooled12 <- p2
    
    # Z-statistics for the two treatment durations
    
    z1 <-  (p_dur11-p_plac11+noninf)/sqrt(sdpooled11)
    
    z2 <- (p_dur12-p_plac11+noninf)/sqrt(sdpooled12)
    
    
    # Count when null hypothesis are rejected at the end of the stage
    
    
    if( z1 >= func$upperbounds[1] &  z2 >= func$upperbounds[1]){
      
      rej_h01andh02[[1]][k,j] = TRUE
      
      tfpc[k,j] = time_interim_month
      
    }
    
    
    if( z1 >= func$upperbounds[1] &  z2 < func$lowerbounds[1]){
      
      rej_h01noh02[[1]][k,j] = TRUE
      
      tfpc[k,j] = time_interim_month
      
      
    }
    
    if( z1 < func$lowerbounds[1] &  z2 >= func$upperbounds[1]){
      
      rej_h02noh01[[1]][k,j] = TRUE
      
      tfpc[k,j] = time_interim_month
      
      
    }
    if( z1 <= func$lowerbounds[1] &  z2 <= func$lowerbounds[1]){
      
      bothineff_all[k,j] = TRUE
      tfpc[k,j] = time_interim_month#+days_over + durations[1]*30.41
      
      
    }
    
    # Continue with both treatment arms
    
    if(z1 < func$upperbounds[1]  & z1 > func$lowerbounds[1]  & z2 > func$lowerbounds[1] &
       z2 < func$upperbounds[1] ){  
      
      plac21 <- plac1st
      dur21 <- dur1st
      dur22 <- dur2st
      
      # continue to recruit until observe the patients on the control duration
      
      plac_2missing <- if_else((2*popord-(length(plac21)))>0, 
                               (2*popord-(length(plac21))), -1)
      
      plac_2 <- months_2 <- 0
      
      count_2 <- NULL
      
      #recruit until get popord on the shortest
      
      while(plac_2 < plac_2missing){
        count_2 <- c(count_2,sample(c("C","L","S"), recrate, prob = p, replace=TRUE))
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
      
      dur22 <- c(dur22, rbinom(n = n_short_2, prob= mean2, size = 1))
      
      #duration
      
      dur[k,j] = dur[k,j]+(months_2 + durations[1])*(months_2>0)+(0)*(months_2==0)#+6
      
      
      
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
      
      
      if (z21 >= func$upperbounds[2] & z22 >= func$upperbounds[2]){
        
        rej_h01andh02[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
        
      }
      
      if (z21 >= func$upperbounds[2] & z22 < func$upperbounds[2]){
        
        rej_h01noh02[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
        
      }
      
      if (z21 < func$upperbounds[2] & z22 >= func$upperbounds[2]){
        
        rej_h02noh01[[2]][k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
        
      }
      if( z21 <= func$upperbounds[2] &  z22 <= func$upperbounds[2]){
        
        bothineff_all[k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
        
      }
      
    }
    
    # Continue with only first treatment arm
    
    if((z1 < func$upperbounds[1] &   z1 > func$lowerbounds[1] & z2 <= func$lowerbounds[1] )|| 
       (z1 < func$upperbounds[1] & z1 > func$lowerbounds[1] & z2 >= func$upperbounds[1])){
      
      plac21 <- plac1st
      dur21 <- dur1st
      
      # continue to recruit until observe the patients on the control duration
      
      plac_2missing <- if_else((2*popord-(length(plac21)))>0, 
                               (2*popord-(length(plac21))), -1)
      
      plac_2 <- months_2 <- 0
      
      count_2 <- NULL
      
      #recruit until get popord on the shortest
      
      while(plac_2 < plac_2missing){
        count_2 <- c(count_2,sample(c("C","L"), recrate, prob = p[c(1,2)], replace=TRUE))
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
      
      estimatedss[k,j] <- length(plac21)+length(dur21)+length(dur2st)
      
      
      # Estimated response rates
      
      p_plac21 <- mean(plac21)
      
      p_dur21 <- mean(dur21)
      
      p21 <- (mean1*(1-mean1)/length(dur21))+(mean0*(1-mean0)/length(plac21))
      
      sdpooled21 <- p21
      
      # Z-statistics for the two treatment durations
      
      z21 <- (p_dur21-p_plac21+noninf)/sqrt(sdpooled21)
      
      # Rejection of the null hypothesis at the second stage
      
      if (z21 >= func$upperbounds[2] & z2 >= func$upperbounds[1]){
        
        rej_h01andh02[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
      }
      
      if (z21 >= func$upperbounds[2] & z2 <= func$lowerbounds[1]){
        
        rej_h01noh02[[2]][k,j] = TRUE
        
        tfpc[k,j] = dur[k,j]
        
      }
      
      if( z21 <= func$upperbounds[2] & z2 >= func$upperbounds[1]){
        
        rej_h02noh01[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
        
      }
      
      if( z21 <= func$upperbounds[2] & z2 <= func$lowerbounds[1]){
        
        bothineff_all[k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
        
      }
    }
    
    # Continue with only second treatment arm
    
    if((z2 < func$upperbounds[1] & z2 > func$lowerbounds[1] & z1 <= func$lowerbounds[1])|| 
       (z2 < func$upperbounds[1] & z2 > func$lowerbounds[1] & z1 >= func$upperbounds[1])){ 
      
      plac21 <- plac1st
      
      dur22 <- dur2st
      
      # continue to recruit until observe the patients on the control duration
      
      plac_2missing <- if_else((2*popord-(length(plac21)))>0, 
                               (2*popord-(length(plac21))), -1)
      
      plac_2 <- months_2 <- 0
      
      count_2 <- NULL
      
      #recruit until get popord on the shortest
      
      while(plac_2 < plac_2missing){
        count_2 <- c(count_2,sample(c("C","S"), recrate, prob = p[c(1,3)], replace=TRUE))
        months_2 <- length(count_2)/recrate
        plac_2 <- sum(count_2=="C")
      }
      #table(count_2)
      
      n_plac_2 <- sum(count_2=="C")
      
      n_short_2 <- sum(count_2=="S")
      
      # Simulate placebo and treatment populations
      
      plac21 <- c(plac21, rbinom(n = n_plac_2, prob= mean0, size = 1))
      
      
      dur22 <- c(dur22, rbinom(n = n_short_2, prob= mean2, size = 1))
      
      #duration
      
      dur[k,j] = dur[k,j]+(months_2 + durations[1])*(months_2>0)+(0)*(months_2==0)#+6
      
      #ESS
      
      estimatedss[k,j] <- length(plac21)+length(dur22)+length(dur1st)
      
      
      # Estimated response rates
      
      p_plac21 <- mean(plac21)
      
      p_dur22 <- mean(dur22)
      
      p22 <- (mean2*(1-mean2)/length(dur22))+(mean0*(1-mean0)/length(plac21))
      
      
      sdpooled22 <- p22
      
      # Z-statistics for the two treatment durations
      
      z22 <- (p_dur22-p_plac21+noninf)/sqrt(sdpooled22)
      
      # Rejection of the null hypothesis at the second stage
      
      if (z22 >= func$upperbounds[2] & z1 >= func$upperbounds[1]){
        
        
        rej_h01andh02[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
      }
      
      if (z22 >= func$upperbounds[2] & z1 <= func$lowerbounds[1]){
        
        rej_h02noh01[[2]][k,j] = TRUE
        
        tfpc[k,j] = dur[k,j]
        
      }
      
      if( z22 <= func$upperbounds[2] & z1 >= func$upperbounds[1]){
        
        rej_h01noh02[[2]][k,j] = TRUE
        tfpc[k,j] = time_interim_month
        
        
      }
      
      if( z22 <= func$upperbounds[2] & z1 <= func$lowerbounds[1]){
        
        bothineff_all[k,j] = TRUE
        tfpc[k,j] = dur[k,j]
        
        
      }
      
      
      
    }
  }
  
  
  # Average Probability of rejecting the null hypotheses
  
  prop_rej[[1]][k,1] <- sum(rej_h01andh02[[1]][k,])/nsim
  
  prop_rej[[1]][k,2] <- sum(rej_h01noh02[[1]][k,])/nsim
  
  prop_rej[[1]][k,3] <- sum(rej_h02noh01[[1]][k,])/nsim
  
  prop_rej[[2]][k,1] <- sum(rej_h01andh02[[2]][k,])/nsim
  
  prop_rej[[2]][k,2] <- sum(rej_h01noh02[[2]][k,])/nsim
  
  prop_rej[[2]][k,3] <- sum(rej_h02noh01[[2]][k,])/nsim
  
  prop_rej[[2]][k,4] <- sum(bothineff_all[k,])/nsim
  
  estimatedsamplesize[k] <- mean(estimatedss[k,]) 
  
  durationtrial[k] <- mean(dur[k,])
  
  timetofirstposclaim[k] <- mean(tfpc[k,])
}

colnames(prop_rej[[1]]) <- c(
  "Allprom_rejH01andH02_1stage", 
  "Allprom_rejH01noH02_1stage",
  "Allprom_rejH02noH01_1stage")

colnames(prop_rej[[2]]) <- c(
  "Allprom_rejH01andH02_2stage", 
  "Allprom_rejH01noH02_2stage",
  "Allprom_rejH02noH01_2stage",
  "Allprom_rejneitherH01andH02")

colnames(estimatedsamplesize) <- "ESS All promising"

colnames(durationtrial) <- "Expected duration (months)"

colnames(timetofirstposclaim) <- "Expected TFPC (months)"

summary_MAMS <- cbind(popord,
                      sum(ceiling((rho*2*popord))),
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
  
  
  
  Allprom_rejH01noH02 = Allprom_rejH01noH02_1stage+Allprom_rejH01noH02_2stage,
  
  Allprom_rejH02noH01 = Allprom_rejH02noH01_1stage+Allprom_rejH02noH01_2stage,
  
  Allprom_rejH01andH02 = Allprom_rejH01andH02_1stage+Allprom_rejH01andH02_2stage,
  
  sumprob = Allprom_rejH01noH02+Allprom_rejH02noH01+Allprom_rejH01andH02+Allprom_rejneitherH01andH02,
  
  typeI = Allprom_rejH01noH02+Allprom_rejH02noH01+Allprom_rejH01andH02,
  
  `reject all` = Allprom_rejH01andH02,
  
  `reject at least one` = typeI,
  
  u1 = func$upperbounds[1],
  u2 = func$upperbounds[2],
  l1 = func$lowerbounds[1]
  
) %>% mutate(
  
  prob = paste("(", p0, ", ", pL,", ",pS, ")", sep = ""),
  theta = paste("(", round(thetaL, 3), ", ", round(thetaS, 3), ")", sep = "")
)
# )%>% select(
#   
#   `patientscontrol_perstage ORD`,
#   `Max theoretical SS`,
#   `Actual SS`,
#   `ESS All promising`,
#   prob,
#   theta,
#   typeI,
#   `reject all`,
#   `reject at least one`,
#   `Expected duration (days)`,
#   `Expected duration (months)`,
#   `Expected TFPC (days)`,
#   `Expected TFPC (months)`,
#   sumprob,
#   u1,
#   u2,
#   l1
# )

write.csv(summary_MAMS, paste("3arm2stageMAMS_TFPC_new3_recrutrate",  recrate, "_", Sys.Date(), 
                              ".csv", sep = ""),
           row.names = F) 

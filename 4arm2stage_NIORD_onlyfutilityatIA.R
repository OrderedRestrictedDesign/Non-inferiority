rm(list=ls())

##### Time to first positive claim

# Durations of the treatment

durations <- c(6,4,3,2)

# Recruitment rate per month

recrate <- 30

############# calls functions (4-arm 2-stage design) with binary data

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

bounds_4arm2stageORDSI <- function(u, arms,ratio,
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
    
    
    prob <-  pmvnorm(lower = c(-Inf,-Inf,l1, u2,-Inf, -Inf),
                     upper = c(Inf,Inf,Inf, Inf,Inf,Inf),
                     sigma = cov)[1]+
      pmvnorm(lower = c(-Inf,l1,-Inf, u2,-Inf, -Inf),
              upper = c(Inf,Inf,l1,Inf,Inf,Inf),
              sigma = cov)[1]+
      pmvnorm(lower = c(l1,-Inf,-Inf, u2,-Inf, -Inf),
              upper = c(Inf,l1,l1,Inf, Inf, Inf),
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

boundaries_4arm2stageORDNonInf_Zmod_SI <- function(theta,
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
  
  rhonull1 <- c(1,1,1,1)
  
  # 4-arm covariance matrix
  
  covbfin12 <- (sqrt(rhonull1[1]*rhonull1[2])*pc*(1-pc))/(sqrt(rhonull1[4]*vart1+rhonull1[1]*varpc)*sqrt(rhonull1[4]*vart2+rhonull1[2]*varpc))
  
  covbfin13 <- (sqrt(rhonull1[1]*rhonull1[3])*pc*(1-pc))/(sqrt(rhonull1[4]*vart1+rhonull1[1]*varpc)*sqrt(rhonull1[4]*vart3+rhonull1[3]*varpc))
  
  covbfin14 <- sqrt(rhonull1[1]*rhonull1[4]*r[2]*rhonull[1]*r[2]*rhonull[4])/sqrt((r[2]*rhonull[4]*vart1+r[2]*rhonull[1]*varpc)*(r[1]*rhonull1[4]*vart1+r[1]*rhonull1[1]*varpc))*((vart1/(r[2]*rhonull[1]))+(varpc/(r[2]*rhonull[4])))
  
  covbfin15 <- sqrt(rhonull1[1]*rhonull1[4]*r[2]*rhonull[2]/(r[2]*rhonull[4]))*varpc/(sqrt(rhonull1[4]*vart1+rhonull1[1]*varpc)*sqrt(r[2]*rhonull[4]*vart2+r[2]*rhonull[2]*varpc))
  
  covbfin16 <- sqrt(rhonull1[1]*rhonull1[4]*r[2]*rhonull[3]/(r[2]*rhonull[4]))*varpc/(sqrt(rhonull1[4]*vart1+rhonull1[1]*varpc)*sqrt(r[2]*rhonull[4]*vart3+r[2]*rhonull[3]*varpc))
  
  covbfin23 <- (sqrt(rhonull1[2]*rhonull1[3])*pc*(1-pc))/(sqrt(rhonull1[4]*vart2+rhonull1[2]*varpc)*sqrt(rhonull1[4]*vart3+rhonull1[3]*varpc))
  
  covbfin24 <- sqrt(rhonull1[2]*rhonull1[4]*r[2]*rhonull[1]/(r[2]*rhonull[4]))*varpc/(sqrt(rhonull1[4]*vart2+rhonull1[2]*varpc)*sqrt(r[2]*rhonull[4]*vart1+r[2]*rhonull[1]*varpc))
  
  covbfin25 <- sqrt(rhonull1[2]*rhonull1[4]*r[2]*rhonull[2]*r[2]*rhonull[4])/sqrt((r[2]*rhonull[4]*vart2+r[2]*rhonull[2]*varpc)*(r[1]*rhonull1[4]*vart2+r[1]*rhonull1[2]*varpc))*((vart2/(r[2]*rhonull[2]))+(varpc/(r[2]*rhonull[4])))
  
  covbfin26 <- sqrt(rhonull1[2]*rhonull1[4]*r[2]*rhonull[3]/(r[2]*rhonull[4]))*varpc/(sqrt(rhonull1[4]*vart2+rhonull1[2]*varpc)*sqrt(r[2]*rhonull[4]*vart3+r[2]*rhonull[3]*varpc))
  
  covbfin34 <- sqrt(rhonull1[3]*rhonull1[4]*r[2]*rhonull[1]/(r[2]*rhonull[4]))*varpc/(sqrt(rhonull1[4]*vart3+rhonull1[3]*varpc)*sqrt(r[2]*rhonull[4]*vart1+r[2]*rhonull[1]*varpc))
  
  covbfin35 <- sqrt(rhonull1[3]*rhonull1[4]*r[2]*rhonull[2]/(r[2]*rhonull[4]))*varpc/(sqrt(rhonull1[4]*vart3+rhonull1[3]*varpc)*sqrt(r[2]*rhonull[4]*vart2+r[2]*rhonull[2]*varpc))
  
  covbfin36 <- sqrt(rhonull1[3]*rhonull1[4]*r[2]*rhonull[3]*r[2]*rhonull[4])/sqrt((r[2]*rhonull[4]*vart3+r[2]*rhonull[3]*varpc)*(r[1]*rhonull1[4]*vart3+r[1]*rhonull1[3]*varpc))*((vart3/(r[2]*rhonull[3]))+(varpc/(r[2]*rhonull[4])))
  
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
  
  first <- bounds_4arm2stageORDSI(u11, arms,r, 
                                  alpha, covbfin, firststagematrix,
                                  ushape, lshape,prec)
  
  low <- first$alow
  
  up <- first$aup
  
  p <- rep(1, times = prec)
  
  for (p in 1:length(p)){
    
    callf <- bounds_4arm2stageORDSI(seq(from = low, 
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
  
  
  covbfin12 <- (sqrt(rhonull1[1]*rhonull1[2])*pc*(1-pc))/(sqrt(rhonull1[4]*vart1+rhonull1[1]*varpc)*sqrt(rhonull1[4]*vart2+rhonull1[2]*varpc))
  
  covbfin13 <- (sqrt(rhonull1[1]*rhonull1[3])*pc*(1-pc))/(sqrt(rhonull1[4]*vart1+rhonull1[1]*varpc)*sqrt(rhonull1[4]*vart3+rhonull1[3]*varpc))
  
  covbfin14 <- sqrt(rhonull1[1]*rhonull1[4]*r[2]*rhopower[1]*r[2]*rhopower[4])/sqrt((r[2]*rhopower[4]*vart1+r[2]*rhopower[1]*varpc)*(r[1]*rhonull1[4]*vart1+r[1]*rhonull1[1]*varpc))*((vart1/(r[2]*rhopower[1]))+(varpc/(r[2]*rhopower[4])))
  
  covbfin15 <- sqrt(rhonull1[1]*rhonull1[4]*r[2]*rhopower[2]/(r[2]*rhopower[4]))*varpc/(sqrt(rhonull1[4]*vart1+rhonull1[1]*varpc)*sqrt(r[2]*rhopower[4]*vart2+r[2]*rhopower[2]*varpc))
  
  covbfin16 <- sqrt(rhonull1[1]*rhonull1[4]*r[2]*rhopower[3]/(r[2]*rhopower[4]))*varpc/(sqrt(rhonull1[4]*vart1+rhonull1[1]*varpc)*sqrt(r[2]*rhopower[4]*vart3+r[2]*rhopower[3]*varpc))
  
  covbfin23 <- (sqrt(rhonull1[2]*rhonull1[3])*pc*(1-pc))/(sqrt(rhonull1[4]*vart2+rhonull1[2]*varpc)*sqrt(rhonull1[4]*vart3+rhonull1[3]*varpc))
  
  covbfin24 <- sqrt(rhonull1[2]*rhonull1[4]*r[2]*rhopower[1]/(r[2]*rhopower[4]))*varpc/(sqrt(rhonull1[4]*vart2+rhonull1[2]*varpc)*sqrt(r[2]*rhopower[4]*vart1+r[2]*rhopower[1]*varpc))
  
  covbfin25 <- sqrt(rhonull1[2]*rhonull1[4]*r[2]*rhopower[2]*r[2]*rhopower[4])/sqrt((r[2]*rhopower[4]*vart2+r[2]*rhopower[2]*varpc)*(r[1]*rhonull1[4]*vart2+r[1]*rhonull1[2]*varpc))*((vart2/(r[2]*rhopower[2]))+(varpc/(r[2]*rhopower[4])))
  
  covbfin26 <- sqrt(rhonull1[2]*rhonull1[4]*r[2]*rhopower[3]/(r[2]*rhopower[4]))*varpc/(sqrt(rhonull1[4]*vart2+rhonull1[2]*varpc)*sqrt(r[2]*rhopower[4]*vart3+r[2]*rhopower[3]*varpc))
  
  covbfin34 <- sqrt(rhonull1[3]*rhonull1[4]*r[2]*rhopower[1]/(r[2]*rhopower[4]))*varpc/(sqrt(rhonull1[4]*vart3+rhonull1[3]*varpc)*sqrt(r[2]*rhopower[4]*vart1+r[2]*rhopower[1]*varpc))
  
  covbfin35 <- sqrt(rhonull1[3]*rhonull1[4]*r[2]*rhopower[2]/(r[2]*rhopower[4]))*varpc/(sqrt(rhonull1[4]*vart3+rhonull1[3]*varpc)*sqrt(r[2]*rhopower[4]*vart2+r[2]*rhopower[2]*varpc))
  
  covbfin36 <- sqrt(rhonull1[3]*rhonull1[4]*r[2]*rhopower[3]*r[2]*rhopower[4])/sqrt((r[2]*rhopower[4]*vart3+r[2]*rhopower[3]*varpc)*(r[1]*rhonull1[4]*vart3+r[1]*rhonull1[3]*varpc))*((vart3/(r[2]*rhopower[3]))+(varpc/(r[2]*rhopower[4])))
  
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
    
    sigma1 <- (varpc/(rhonull1[4]*pop[n]))+(vart1/(rhonull1[1]*pop[n]))
    
    sigma2 <- (varpc/(rhonull1[4]*pop[n]))+(vart2/(rhonull1[2]*pop[n]))
    
    sigma3 <- (varpc/(rhonull1[4]*pop[n]))+(vart3/(rhonull1[3]*pop[n]))
    
    sigma21 <- (varpc/(rhopower[4]*r[2]*pop[n]))+(vart1/(rhopower[1]*r[2]*pop[n]))
    
    sigma22 <- (varpc/(rhopower[4]*r[2]*pop[n]))+(vart2/(rhopower[2]*r[2]*pop[n]))
    
    sigma23 <- (varpc/(rhopower[4]*r[2]*pop[n]))+(vart3/(rhopower[3]*r[2]*pop[n]))
    
    sigman <- c(sigma1,sigma2,sigma3,sigma21,sigma22,sigma23)
    
    mean1 <- rep(((theta[2:(arms)]-theta[1])+noninf),2)/(sqrt(sigman))
    
    if(power == "Reject all hyp."){
      
      set.seed(seed)
      
      prob <- pmvnorm(lower = c(-Inf,-Inf,l1, u2, u2, u2),
                      upper = c(Inf, Inf, Inf, Inf, Inf, Inf),
                      sigma = covbfinB,
                      mean = mean1)[1]
    }
    if(power == "Reject at least one hyp."){
      
      set.seed(seed)
      
      prob <- pmvnorm(lower = c(-Inf,-Inf,l1, u2,-Inf, -Inf),
                      upper = c(Inf,Inf,Inf, Inf,Inf,Inf),
                      sigma = covbfinB,
                      mean = mean1)[1]+
        pmvnorm(lower = c(-Inf,l1,-Inf, u2,-Inf, -Inf),
                upper = c(Inf,Inf,l1,Inf,Inf,Inf),
                sigma = covbfinB,
                mean = mean1)[1]+
        pmvnorm(lower = c(l1,-Inf,-Inf, u2,-Inf, -Inf),
                upper = c(Inf,l1,l1,Inf, Inf, Inf),
                sigma = covbfinB,
                mean = mean1)[1]
    }
    
    if(power == "Reject LongANDMedium"){
      
      set.seed(seed)
      
      prob <- pmvnorm(lower = c(-Inf,-Inf,l1, u2, u2, -Inf),
                      upper = c(Inf, Inf, Inf, Inf, Inf, Inf),
                      sigma = covbfinB,
                      mean = mean1)[1]+
        pmvnorm(lower = c(-Inf,l1,-Inf, u2, u2, -Inf),
                upper = c(Inf, Inf, l1, Inf, Inf, Inf),
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

#### Simulations #####


# PARAMETERS for both designs

seed <- 64736 # set seed to simulate the data

nsim <- 10^5 # number of simulations

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

r <-  c(1,2)

# call function to find the sample size per arm per stage and the critical boundaries

fun <- boundaries_4arm2stageORDNonInf_Zmod_SI(theta = c(pc,pt1,pt2,pt3),
                                              noninf = noninf,
                                              stage = stage,
                                              alpha = alpha,
                                              beta = beta,
                                              r = r,
                                              rhonull = rhopower,
                                              rhopower = rhopower,
                                              prec = prec,
                                              arms = arms,
                                              ushape = ushape,
                                              lshape = ushape,
                                              seed = seed,
                                              power=powerreq)
fun

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

scenario <- scenario[c(8,9),]

######## Sample size and bound grid values 

popord <- fun$`sample size per arm per stage`


# triangular bounds

summary_ORD <- NULL

u2 <- fun$upperbounds1[2]

l1 <- fun$lowerbounds1


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

rej_ord_h01noh02noh03 <- matrix(data=FALSE, nrow=nrow(scenario),
                                ncol=nsim)

rej_ord_h01andh02noh03 <- matrix(data=FALSE, nrow=nrow(scenario),
                                 ncol=nsim)

prop_ord <- list()

prop_ord[[1]] <- matrix(data=NA, ncol=2,
                        nrow = nrow(scenario))
prop_ord[[2]] <- matrix(data=NA, ncol=4,
                        nrow = nrow(scenario))

estimatedss <- actualmaxss <- matrix(data=NA, nrow=nrow(scenario),
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
      count <- c(count,sample(c("C","L","M", "S"), recrate, prob = p, replace=TRUE))
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
    
    time_interim_month <- max(c(month_lastpatplac+durations[1],
                                month_lastpatlong+durations[2],
                                month_lastpatmedium+durations[3],
                                month_lastpatshort+durations[4]))
    
    index_timeinter <- which.max(c(month_lastpatplac+durations[1],
                                   month_lastpatlong+durations[2],
                                   month_lastpatmedium+durations[3],
                                   month_lastpatshort+durations[4]))
    
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
    
    
    if(z1 <= l1 &  z2 <= l1  
       & z3 <= l1){
      
      rej_any[[1]][k,j] = TRUE
      
      
    }
    
    # Continue to second stage with all arms
    
    if(( z3 > l1)){
      
      plac21 <- plac1st
      dur21 <- dur1st
      dur22 <- dur2st
      dur23 <- dur3st
      
      # continue to recruit until observe the patients on the control duration
      
      plac_2missing <-  if_else((2*popord-(length(plac21)))>0, 
                                (2*popord-(length(plac21))), -1)
      
      plac_2 <- months_2 <-0
      
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
      
      dur[k,j] = dur[k,j]+(months_2 + durations[1])*(months_2>0)+(0)*(months_2==0)
      
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
      
      
      if( z21 >= u2 &  
          z22 >= u2 &
          z23 >= u2){
        
        rej_all[[2]][k,j] = TRUE
        
      }
      
      if( z21 < u2){
        
        rej_any[[2]][k,j] = TRUE
        
        
      }
      
      if(  z21 >= u2 &  
           z22 >= u2 &
           z23 < u2){
        
        rej_ord_h01andh02noh03[k,j] = TRUE
        
      }
      
      if( z21 >= u2 &  
          z22 < u2){
        
        rej_ord_h01noh02noh03[k,j] = TRUE
        
      }
      
    }
    # Continue to second stage with arms 1 and 2
    
    if(
      ( z2 > l1  &  z3 < l1)
    ){
      
      # Simulate placebo and treatment populations
      
      plac21 <- plac1st
      dur21 <- dur1st
      dur22 <- dur2st
      
      # continue to recruit until observe the patients on the control duration
      
      plac_2missing <-if_else((2*popord-(length(plac21)))>0, 
                              (2*popord-(length(plac21))), -1)
      
      plac_2 <- months_2 <-0
      
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
      
      dur[k,j] =  dur[k,j]+(months_2 + durations[1])*(months_2>0)+(0)*(months_2==0)
      
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
      
      # Rejection of the null hypothesis at the second stage
      
      if (z21 < u2){
        
        rej_any[[2]][k,j]  = TRUE
        
      }
      
      if(  z21 >= u2 &  
           z22 >= u2){
        
        rej_ord_h01andh02noh03[k,j] = TRUE
        
      }
      
      if(  z21 >= u2 &  
           z22 < u2  ){
        
        rej_ord_h01noh02noh03[k,j] = TRUE
        
      }
      
    }
    
    
    # Continue to second stage with 1
    
    if((z1 > l1  & z2 <  l1 &  z3 < l1)
    ){
      
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
      
      dur[k,j] = dur[k,j]+(months_2 + durations[1])*(months_2>0)+(0)*(months_2==0)
      #ESS
      
      
      estimatedss[k,j] <- length(plac21)+length(dur21)+length(dur2st)+length(dur3st)
      
      
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
  
  durationtrial[k] <- mean(dur[k,])
  
}



colnames(prop_ord[[1]]) <- c("ORD_rejall_1stage",
                             "ORD_rejany_1stage")

colnames(prop_ord[[2]]) <- c("ORD_rejall_2stage",
                             "ORD_rejany_2stage",
                             "ORD_rejH01andH02noH03",
                             "ORD_rejH01noH02noH03")

colnames(estimatedsamplesize) <- "ESS Ordered ORD"

colnames(durationtrial) <- "Expected duration (months)"

summary_ORD <- cbind(popord,
                     sum(ceiling((rhopower*2*popord))),
                     apply(actualmaxss,1,mean, na.rm=TRUE),
                     scenario,
                     prop_ord,
                     estimatedsamplesize,
                     durationtrial)

colnames(summary_ORD)[1:3]<- c("patientscontrol_perstage ORD",
                               "Max theoretical SS",
                               "Actual SS")

summary_ORD <- summary_ORD %>% mutate(
  
  rej_all_overall = ORD_rejall_1stage + ORD_rejall_2stage,
  
  rej_any_overall = ORD_rejany_1stage + ORD_rejany_2stage,
  
  rej_LongANDMedium = ORD_rejH01andH02noH03 + rej_all_overall,
  
  rej_LongOnly = ORD_rejH01noH02noH03, 
  
  typeI = ORD_rejH01noH02noH03 + ORD_rejH01andH02noH03 + rej_all_overall,
  
  sumprob = typeI + rej_any_overall
  
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
  typeI,
  rej_LongANDMedium,
  rej_all_overall,
  rej_any_overall,
  `Expected duration (months)`
  
  
)

summary_ORD <- cbind(summary_ORD, u2,l1)

write.csv(summary_ORD, paste("4arm2stage_onlyfutility_recrate",recrate, "_", Sys.Date(),
".csv", sep = ""),row.names = F)

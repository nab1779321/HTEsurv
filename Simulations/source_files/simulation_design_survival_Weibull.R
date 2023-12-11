expit <- function(x) {
  exp(x) / (1 + exp(x))
}

rm(list=ls())

data_gen_censor = function(n=5000, p=10,  PH=TRUE, censor = "30%", setting = 1, ex = "random",time.interest=NULL) {
  X <- matrix(rnorm(p * n, sd = 0.35), nrow = n, ncol = p)
  x1 <- X[, 1]
  x2 <- X[, 2]
  x3 <- X[, 3]
  x4 <- X[, 4]
  x5 <- X[, 5]
  x6 <- X[, 6]
  x7 <- X[, 7]
  x8 <- ifelse(X[, 8]<=0,0,1)
  x9 <- ifelse(X[, 9]<=0,0,1)
  x10 <- ifelse(X[, 10]<=0,0,1)
  
  if (ex == "random") {
    p1 <- rep(0.5,n)
  }
  
  
  if (ex == "unbalanced") {
    p1 <- rep(0.05,n)
  }
  
  if (ex == "dependent") {
    p1 <- expit(
      1.3*x1-0.8*x5
    )
  }
  
  
  p0 <- 1 - p1
  Trt = NULL
  for (i in 1:n) {
    Trt[i] <- sample(c(1, 0),
                     size = 1,
                     replace = TRUE,
                     prob = c(p1[i], p0[i]))
  }
  
  if (PH == FALSE) {
    eta_1 <- exp(0.7 + 1.8 * x3 + 0.8 * x7)
    eta_2 <- exp(0.9 - 0.5 * x1 + 0.5 * x2)
    
    summary(eta_1)
  } 
  if (PH == TRUE) {
    eta_1 <- 2
    eta_2 <- 2
  }
  
  if (setting == 1) {
    ##median time = 15
    LP1 <- 0.2 * x1 + 0.7 * x2  + 0.4 * x9
    LP2 <- - 0.5 * x1 - 2 * x2  - 0.25 * x9
  }
  
  if (setting == 2) {
    ##median time = 15
    LP1 <- -0.5*x1+0.7*x2+0.2*x9+0.9*x2*x9
    LP2 <- -0.05* exp(x1)-0.2*x2^2+0.35*x9
  }
  
  if (setting == 3) {
    ##median time = 14
    LP1 <- -0.5*x1+0.7*x2+0.2*x9+0.9*x2*x9+0.6*x3-0.5*x4^2+0.6*x8
    LP2 <- -0.05*exp(x1)-0.2*x2^2 +0.2*x9-0.1*exp(x5)+0.7*sin(x6)+0.5*x10
  }
  
  #independent censoring
  if (censor == "30%" & setting==1){
    C <- rexp(n, rate = 0.02)
  }
  
  if (censor == "30%" & setting!=1){
    C <- rexp(n, rate = 0.025)
  }
  
  #generate U ~ unif(0,1)
  U = runif(n, 0,1)
  #scale parameter lambda>0, exp(Linear predictor)
  lambda1 <- 18
  lambda2 <- 20
  
  #potential survival times
  T1 <- lambda1*(((-log(U))/exp(LP1))^(1/eta_1))
  T2 <- lambda2*(((-log(U))/exp(LP2))^(1/eta_2))
  # observed outcomes
  T <- cbind(T1, T2)
  TTrt <- cbind(T, Trt)
  Tobs <- apply(TTrt, 1, function(x)
    x[1:2][x[3] + 1]) #T2 is observed if Trt=1, and T1 is observed if Trt=0
  
  
  Tobs_C <- pmin(Tobs, C)
  
  #censoring indicator
  delta <- ifelse(Tobs > C, 0, 1)
  #summary(Tobs_C)
  
  dat<-data.frame(cbind(V1=x1,V2=x2,V3=x3,V4=x4,V5=x5,V6=x6,V7=x7,V8=x8,V9=x9,V10=x10,Treatment=Trt,Time = Tobs_C, Event=delta))
  
  # True survival probability at time.interest
  time.interest<-time.interest
  true_mst_trt_1 <- exp(-exp(LP1) * ((time.interest/lambda1) ^ eta_1))
  true_mst_trt_2 <- exp(-exp(LP2) * ((time.interest/lambda2) ^ eta_2))
  true.diff<-true_mst_trt_2-true_mst_trt_1
  true.ratio<-true_mst_trt_2/true_mst_trt_1
  
  
  return(
    list(
      data = dat,
      true.diff = true.diff,
      true.ratio = true.ratio,
      true1=true_mst_trt_2,
      true0=true_mst_trt_1
    )
  )
  
}

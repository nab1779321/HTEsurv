expit <- function(x) {
  exp(x) / (1 + exp(x))
}

data_gen_censor_slogistic= function(n=1000, p=10,sigma0=1,sigma1=2,r=0.15,mu=0, censor = "30%", setting = 1, ex = "random",time.interest=1.05) { 
  
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
  
    if (setting == 1) {
      ##median time = 1.05
      LP0 <- 0.2 * x1 + 0.7 * x2  + 0.4 * x9
      LP1 <- - 0.5 * x1 - 2 * x2  - 0.25 * x9 
    }
    
    if (setting == 2) { 
      ##median time = 1.1
      LP0 <- -0.5*x1+0.7*x2+0.2*x9+0.9*x2*x9
      LP1 <- -0.05* exp(x1)-0.2*x2^2+0.35*x9
    }
    
    if (setting == 3) {
      ##median time = 1.29
      LP0 <- -0.5*x1+0.7*x2+0.2*x9+0.9*x2*x9+0.6*x3-0.5*x4^2+0.6*x8
      LP1 <- -0.05*exp(x1)-0.2*x2^2 +0.2*x9-0.1*exp(x5)+0.7*sin(x6)+0.5*x10
    }
  
  #independent censoring
  if (censor == "30%" & setting==1){
    C <- rexp(n, rate = r)  # r=0.15
  }
  
  
  if (censor == "30%" & setting==2){ 
    C <- rexp(n, rate = r)  # r=0.22
  }
  
  if (censor == "30%" & setting==3){
    C <- rexp(n, rate = r) # r=0.2
  }
  
  #generate W~logistic(0,1)
  W <- rlogis(n,location = 0,scale = 1)
  
  #potential survival times
  T1 <- exp(mu+LP0+sigma0*W) 
  T0 <- exp(mu+LP1+sigma1*W) 
  # observed outcomes
  T <- cbind(T0, T1)
  TTrt <- cbind(T, Trt)
  Tobs <- apply(TTrt, 1, function(x)
    x[1:2][x[3] + 1]) #T2 is observed if Trt=1, and T1 is observed if Trt=0
  
  
  Tobs_C <- pmin(Tobs, C)
  
  #censoring indicator
  delta <- ifelse(Tobs > C, 0, 1)
  
  dat<-data.frame(cbind(V1=x1,V2=x2,V3=x3,V4=x4,V5=x5,V6=x6,V7=x7,V8=x8,V9=x9,V10=x10,Treatment=Trt,Time = Tobs_C, Event=delta))
  
  # True survival probability at time.interest
  time.interest<-time.interest
  true_mst_trt_1 <- 1/( 1+exp( (log(time.interest)-mu-LP0) /sigma0) ) 
  true_mst_trt_0 <- 1/( 1+exp( (log(time.interest)-mu-LP1) /sigma1) ) 
  true.diff<-true_mst_trt_1-true_mst_trt_0
  
  
  return(
    list(
      data = dat,
      true.diff = true.diff,
      true1=true_mst_trt_1,
      true0=true_mst_trt_0
    )
  )
  
}


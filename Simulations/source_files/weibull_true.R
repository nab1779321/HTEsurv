weibull_true <- function(data=dat,testdat=testdat,time.interest=NULL,setting=1){
  data<-data[complete.cases(data),]
  data1 <- data[which(data$Treatment==1),-c(which(names(data)=="Treatment"))]
  data0 <- data[which(data$Treatment==0),-c(which(names(data)=="Treatment"))]
  
  # Fit Weibull model
  if (setting==1){
    weib0 <- survreg(Surv(data0$Time, data0$Event)~data0$V1+data0$V2+data0$V9,dist="weibull")
    sigma0<-summary(weib0)$scale
    mu0<-summary(weib0)$coef[1]
    lambda0<-exp(mu0)
    eta0<-1/sigma0
    beta0<-summary(weib0)$coef[-1]*(-1)/sigma0
    
    weib1 <- survreg(Surv(data1$Time, data1$Event)~data1$V1+data1$V2+data1$V9,dist="weibull")
    sigma1<-summary(weib1)$scale
    mu1<-summary(weib1)$coef[1]
    lambda1<-exp(mu1)
    eta1<-1/sigma1
    beta1<-summary(weib1)$coef[-1]*(-1)/sigma1
    # Predict using counterfactual data
    time.interest1<-time.interest0<-time.interest
    predict_weib1 <- exp(-exp(as.matrix(testdat[,c(match(c("V1","V2","V9"),names(testdat)))])%*%beta1)*(time.interest/lambda1)^eta1)
    predict_weib0 <- exp(-exp(as.matrix(testdat[,c(match(c("V1","V2","V9"),names(testdat)))])%*%beta0)*(time.interest/lambda0)^eta0)
  }
  
  if (setting==2){
    data0$int<-data0$V2*data0$V9
    weib0 <- survreg(Surv(data0$Time, data0$Event)~data0$V1+data0$V2+data0$V9+data0$int,dist="weibull")
    sigma0<-summary(weib0)$scale
    mu0<-summary(weib0)$coef[1]
    lambda0<-exp(mu0)
    eta0<-1/sigma0
    beta0<-summary(weib0)$coef[-1]*(-1)/sigma0
    
    data1$sqv2<-data1$V2^2
    data1$expv1<-exp(data1$V1)
    weib1 <- survreg(Surv(data1$Time, data1$Event)~data1$expv1+data1$sqv2+data1$V9,dist="weibull")
    sigma1<-summary(weib1)$scale
    mu1<-summary(weib1)$coef[1]
    lambda1<-exp(mu1)
    eta1<-1/sigma1
    beta1<-summary(weib1)$coef[-1]*(-1)/sigma1
    # Predict using counterfactual data
    time.interest1<-time.interest0<-time.interest
    testdat$sqv2<-testdat$V2^2
    testdat$expv1<-exp(testdat$V1)
    testdat$int<-testdat$V2*testdat$V9
    testdat$intercept<-rep(1,nrow(testdat))
    predict_weib1 <- exp(-exp(as.matrix(testdat[,c(match(c("expv1","sqv2","V9"),names(testdat)))])%*%beta1)*(time.interest/lambda1)^eta1)
    predict_weib0 <- exp(-exp(as.matrix(testdat[,c(match(c("V1","V2","V9","int"),names(testdat)))])%*%beta0)*(time.interest/lambda0)^eta0)
  }
  
  if (setting==3){
    data0$int<-data0$V2*data0$V9
    data0$sqv4<-data0$V4^2
    testdat$int<-testdat$V2*testdat$V9
    testdat$sqv4<-testdat$V4^2
    weib0 <- survreg(Surv(data0$Time, data0$Event)~data0$V1+data0$V2+data0$V9+data0$int+data0$V3+data0$sqv4+data0$V8,dist="weibull")
    sigma0<-summary(weib0)$scale
    mu0<-summary(weib0)$coef[1]
    lambda0<-exp(mu0)
    eta0<-1/sigma0
    beta0<-summary(weib0)$coef[-1]*(-1)/sigma0
    
    data1$sqv2<-data1$V2^2
    testdat$sqv2<-testdat$V2^2
    data1$expv1<-exp(data1$V1)
    testdat$expv1<-exp(testdat$V1)
    data1$expv5<-exp(data1$V5)
    testdat$expv5<-exp(testdat$V5)
    data1$sinv6<-sin(data1$V6)
    testdat$sinv6<-sin(testdat$V6)
    weib1 <- survreg(Surv(data1$Time, data1$Event)~data1$expv1+data1$sqv2+data1$V9+data1$expv5+data1$sinv6+data1$V10,dist="weibull")
    sigma1<-summary(weib1)$scale
    mu1<-summary(weib1)$coef[1]
    lambda1<-exp(mu1)
    eta1<-1/sigma1
    beta1<-summary(weib1)$coef[-1]*(-1)/sigma1
    # Predict using counterfactual data
    time.interest1<-time.interest0<-time.interest
    predict_weib1 <- exp(-exp(as.matrix(testdat[,c(match(c("expv1","sqv2","V9","expv5","sinv6","V10"),names(testdat)))])%*%beta1)*(time.interest/lambda1)^eta1)
    predict_weib0 <- exp(-exp(as.matrix(testdat[,c(match(c("V1","V2","V9","int","V3","sqv4","V8"),names(testdat)))])%*%beta0)*(time.interest/lambda0)^eta0)
  }
  
  weib.diff <- predict_weib1-predict_weib0
  weib.ratio <- predict_weib1/predict_weib0
  return(list(diff=weib.diff,ratio=weib.ratio,pred0=predict_weib0,pred1=predict_weib1)) 
}

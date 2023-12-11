slogistic_true <- function(data=dat,testdat=testdat,time.interest=1.05,setting=1){
  data<-data[complete.cases(data),]
  data1 <- data[which(data$Treatment==1),-c(which(names(data)=="Treatment"))]
  data0 <- data[which(data$Treatment==0),-c(which(names(data)=="Treatment"))]
  
    if (setting==1){
      slogistic0 <- survreg(Surv(data0$Time, data0$Event)~data0$V1+data0$V2+data0$V9,dist="loglogistic")
      sigma0<-summary(slogistic0)$scale
      mu0<-summary(slogistic0)$coef[1]
      alpha0 <- 1/sigma0
      lambda0<-exp(-mu0/sigma0)
      beta0<-summary(slogistic0)$coef[-1]*(-1)/sigma0
      
      slogistic1 <- survreg(Surv(data1$Time, data1$Event)~data1$V1+data1$V2+data1$V9,dist="loglogistic")
      sigma1<-summary(slogistic1)$scale
      mu1<-summary(slogistic1)$coef[1]
      alpha1 <- 1/sigma1
      lambda1<-exp(-mu1/sigma1)
      beta1<-summary(slogistic1)$coef[-1]*(-1)/sigma1
      
      # Predict using counterfactual data
      time.interest1<-time.interest0<-time.interest
      predict_slogistic1 <- 1 / ( 1 + exp(as.matrix(testdat[,c(match(c("V1","V2","V9"),names(testdat)))])%*%beta1)*(time.interest1^alpha1)*lambda1 )
      predict_slogistic0 <- 1 / ( 1 + exp(as.matrix(testdat[,c(match(c("V1","V2","V9"),names(testdat)))])%*%beta0)*(time.interest0^alpha0)*lambda0 )
    }
    
    if (setting==2){ 
      data0$sqv2<-data0$V2^2
      data0$expv1<-exp(data0$V1)
      slogistic0 <- survreg(Surv(data0$Time, data0$Event)~data0$expv1+data0$sqv2+data0$V9,dist="loglogistic")
      sigma0<-summary(slogistic0)$scale
      mu0<-summary(slogistic0)$coef[1]
      alpha0 <- 1/sigma0
      lambda0<-exp(-mu0/sigma0)
      beta0<-summary(slogistic0)$coef[-1]*(-1)/sigma0
      
      data1$int<-data1$V2*data1$V9
      slogistic1 <- survreg(Surv(data1$Time, data1$Event)~data1$V1+data1$V2+data1$V9+data1$int,dist="loglogistic")
      sigma1<-summary(slogistic1)$scale
      mu1<-summary(slogistic1)$coef[1]
      alpha1 <- 1/sigma1
      lambda1<-exp(-mu1/sigma1)
      beta1<-summary(slogistic1)$coef[-1]*(-1)/sigma1
      
      #Predict using counterfactual data
      time.interest1<-time.interest0<-time.interest
      testdat$sqv2<-testdat$V2^2
      testdat$expv1<-exp(testdat$V1)
      testdat$int<-testdat$V2*testdat$V9
      predict_slogistic1 <- 1 / ( 1 + exp(as.matrix(testdat[,c(match(c("V1","V2","V9","int"),names(testdat)))])%*%beta1)*(time.interest1^alpha1)*lambda1 )
      predict_slogistic0 <- 1 / ( 1 + exp(as.matrix(testdat[,c(match(c("expv1","sqv2","V9"),names(testdat)))])%*%beta0)*(time.interest0^alpha0)*lambda0 )
    }
    
    if (setting==3){
      data0$sqv2<-data0$V2^2
      testdat$sqv2<-testdat$V2^2
      data0$expv1<-exp(data0$V1)
      testdat$expv1<-exp(testdat$V1)
      data0$expv5<-exp(data0$V5)
      testdat$expv5<-exp(testdat$V5)
      data0$sinv6<-sin(data0$V6)
      testdat$sinv6<-sin(testdat$V6)
      slogistic0 <- survreg(Surv(data0$Time, data0$Event)~data0$expv1+data0$sqv2+data0$V9+data0$expv5+data0$sinv6+data0$V10,dist="loglogistic")
      sigma0<-summary(slogistic0)$scale
      mu0<-summary(slogistic0)$coef[1]
      alpha0 <- 1/sigma0
      lambda0<-exp(-mu0/sigma0)
      beta0<-summary(slogistic0)$coef[-1]*(-1)/sigma0
      
      data1$int<-data1$V2*data1$V9
      data1$sqv4<-data1$V4^2
      testdat$int<-testdat$V2*testdat$V9
      testdat$sqv4<-testdat$V4^2
      slogistic1 <- survreg(Surv(data1$Time, data1$Event)~data1$V1+data1$V2+data1$V9+data1$int+data1$V3+data1$sqv4+data1$V8,dist="loglogistic")
      sigma1<-summary(slogistic1)$scale
      mu1<-summary(slogistic1)$coef[1]
      alpha1 <- 1/sigma1
      lambda1<-exp(-mu1/sigma1)
      beta1<-summary(slogistic1)$coef[-1]*(-1)/sigma1
      
      # Predict using counterfactual data
      time.interest1<-time.interest0<-time.interest
      predict_slogistic1 <- 1 / ( 1 + exp(as.matrix(testdat[,c(match(c("V1","V2","V9","int","V3","sqv4","V8"),names(testdat)))])%*%beta1)*(time.interest1^alpha1)*lambda1 )
      predict_slogistic0 <- 1 / ( 1 + exp(as.matrix(testdat[,c(match(c("expv1","sqv2","V9","expv5","sinv6","V10"),names(testdat)))])%*%beta0)*(time.interest0^alpha0)*lambda0 )
    }
  
  
  slogistic.diff <- predict_slogistic1-predict_slogistic0
  return(list(diff=slogistic.diff,pred0=predict_slogistic0,pred1=predict_slogistic1,param0=param.slogistic0,param1=param.slogistic1)) 
}
rsf_HTE_T <- function(data=dat,testdat=testdat,time.interest=NULL){
  data<-data[complete.cases(data),]
  data1 <- data[which(data$Treatment==1),-c(which(names(data)=="Treatment"))]
  data0 <- data[which(data$Treatment==0),-c(which(names(data)=="Treatment"))]
  
  # Fit RFSRC model
  rfsrc_data1 <-
    rfsrc(
      Surv(Time, Event) ~ .,
      data = data1
    )
  
  rfsrc_data0 <-
    rfsrc(
      Surv(Time, Event) ~ .,
      data = data0
    )
  # Predict using counterfactual data
  time.interest1<-time.interest0<-time.interest
  predict_rfsrc_median1 <- predictSurvProb(rfsrc_data1, newdata = testdat, times=time.interest)
  predict_rfsrc_median0 <- predictSurvProb(rfsrc_data0, newdata = testdat, times=time.interest)
  
  rfsrc.predict0<-predict_rfsrc_median0
  rfsrc.predict1<-predict_rfsrc_median1
  rfsrc.diff <- predict_rfsrc_median1-predict_rfsrc_median0
  return(list(diff=rfsrc.diff,pred0=rfsrc.predict0,pred1=rfsrc.predict1)) 
}

rsf_HTE_X <- function(data=dat,testdat=mydata$data,ntrees=1000,time.interest=15){
  data<-data[complete.cases(data),]
  data1 <- data[which(data$Treatment==1),-c(which(names(data)=="Treatment"))]
  data0 <- data[which(data$Treatment==0),-c(which(names(data)=="Treatment"))]
  
  #Step1: Fit 2 RFSRC models
  rfsrc_data1 <-
    rfsrc(
      Surv(Time, Event) ~ .,
      data = data1
    )
  
  rfsrc_data0 <-
    rfsrc(
      Surv(Time, Event) ~ .,
      data = data0
    )
  time.interest01<-time.interest11<- time.interest00<-time.interest10<-time.interest

    #Step 2: calculate D1 and D2
    # D1=mu1(1)-mu0(1)
    predict_rfsrc_median01 <- predictSurvProb(rfsrc_data0, newdata = data1, times=time.interest01)
    predict_rfsrc_median11 <- predictSurvProb(rfsrc_data1, newdata = data1, times=time.interest11)
    data1$d1<- predict_rfsrc_median11-predict_rfsrc_median01
    data1$r1<- predict_rfsrc_median11/predict_rfsrc_median01
    
    # D0=mu1(0)-mu0(0)
    predict_rfsrc_median00 <- predictSurvProb(rfsrc_data0, newdata = data0, times=time.interest00)
    predict_rfsrc_median10 <- predictSurvProb(rfsrc_data1, newdata = data0, times=time.interest10)
    data0$d0<- predict_rfsrc_median10-predict_rfsrc_median00
    data0$r0<- predict_rfsrc_median10/predict_rfsrc_median00
    
  #Step 3: fit two regressions to model d0 and d1
  modeld0 <-gbm(d0~., data = data0[!is.na(data0$d0),-which(names(data0)%in%c("Time","Event","r0","sim","status"))], 
                distribution="gaussian",n.trees = ntrees,interaction.depth=6, n.minobsinnode = round(959/100),bag.fraction = 0.8)
  
  modeld1 <-gbm(d1~., data = data1[!is.na(data1$d1),-which(names(data1)%in%c("Time","Event","r1","sim","status"))], 
                distribution="gaussian",n.trees = ntrees,interaction.depth=6, n.minobsinnode = round(959/100),bag.fraction = 0.8)
  
  predictiond0<-predict(modeld0,testdat,n.trees=ntrees)
  predictiond1<-predict(modeld1,testdat,n.trees=ntrees)
  
  ####step 4: fit propensity model of treatment assignment
  ####random forest
  data$Treatment<-as.factor(data$Treatment)
  rf<-randomForest(Treatment~., data=data[,-c(ncol(data),ncol(data)-1)],
                   ntree = ntrees, classwt = c(0.5, 0.5))
  pred<-predict(rf,testdat,type="prob")[,2]
  
  X.diff<-pred*predictiond0+(1-pred)*predictiond1
  
  return(list(diff=X.diff)) 

}




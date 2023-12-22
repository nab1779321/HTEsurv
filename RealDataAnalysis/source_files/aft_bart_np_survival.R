baft_HTE_T <- function(data=dat,testdat=testdat,time.interest=NULL){
  data<-data[complete.cases(data),]
  data1 <- data[which(data$Treatment==1),-c(which(names(data)=="Treatment"))]
  data0 <- data[which(data$Treatment==0),-c(which(names(data)=="Treatment"))]

  # Fit AFTrees model 
  AFTrees1 = AFTrees(
    x.train = data1[,-c(which(names(data1)%in%c("Time","Event")))],
    y.train =  data1$Time,
    status = data1$Event,
    nskip = 100,
    ndpost = 1000,
    ntree=500,
    x.test = testdat[,-c(which(names(testdat)%in%c("Time","Event","Treatment")))],
    nonparametric = T
  )
  
  AFTrees0 = AFTrees(
    x.train = data0[,-c(which(names(data0)%in%c("Time","Event")))],
    y.train =  data0$Time,
    status = data0$Event,
    nskip = 100,
    ndpost = 1000,
    ntree=500,
    x.test = testdat[,-c(which(names(testdat)%in%c("Time","Event","Treatment")))],
    nonparametric = T
  )

  # Calculate survival probability
  time.interest1<-time.interest0<-time.interest
  AFTrees1_all_survival_prob <- AFTrees_SurvivalProb (object = AFTrees1, test.only = T,time.points=rep(time.interest1,2))
  AFTrees0_all_survival_prob <- AFTrees_SurvivalProb (object = AFTrees0, test.only = T,time.points=rep(time.interest0,2)) 
  AFTrees1_survival_prob <- AFTrees1_all_survival_prob$Surv.test[,1]
  AFTrees0_survival_prob <- AFTrees0_all_survival_prob$Surv.test[,1]
  

  pred0<-AFTrees0_survival_prob
  pred1<-AFTrees1_survival_prob
  
  baft.diff <- pred1-pred0
  return(list(diff=baft.diff,pred0=pred0,pred1=pred1))
  
}


baft_HTE_X <- function(data=dat,testdat=testdat,ntrees=1000,time.interest=NULL){
  data<-data[complete.cases(data),]
  data1 <- data[which(data$Treatment==1),-c(which(names(data)=="Treatment"))]
  data0 <- data[which(data$Treatment==0),-c(which(names(data)=="Treatment"))]
  
  # Fit AFTrees model 
  AFTrees1 = AFTrees(
    x.train = data1[,-c(which(names(data1)%in%c("Time","Event")))],
    y.train =  data1$Time,
    status = data1$Event,
    nskip = 100,
    ndpost = 1000,
    ntree=500,
    x.test = data[,-c(which(names(data)%in%c("Time","Event","Treatment")))],
    nonparametric = T
  )
  
  AFTrees0 = AFTrees(
    x.train = data0[,-c(which(names(data0)%in%c("Time","Event")))],
    y.train =  data0$Time,
    status = data0$Event,
    nskip = 100,
    ndpost = 1000,
    ntree=500,
    x.test = data[,-c(which(names(data)%in%c("Time","Event","Treatment")))],
    nonparametric = T
  )
  
  
  # Calculate survival probability
  time.interest1<-time.interest0<-time.interest
  AFTrees1_all_survival_prob <- AFTrees_SurvivalProb (object = AFTrees1, test.only = T,time.points=rep(time.interest1,2))
  AFTrees0_all_survival_prob <- AFTrees_SurvivalProb (object = AFTrees0, test.only = T,time.points=rep(time.interest0,2)) 
  AFTrees1_survival_prob <- AFTrees1_all_survival_prob$Surv.test[,2]
  AFTrees0_survival_prob <- AFTrees0_all_survival_prob$Surv.test[,2]
  

  pred0<-AFTrees0_survival_prob
  pred1<-AFTrees1_survival_prob
  
  pred00<-pred0[which(data$Treatment==0)]
  pred01<-pred0[which(data$Treatment==1)]
  pred10<-pred1[which(data$Treatment==0)]
  pred11<-pred1[which(data$Treatment==1)]
  
  data0$d0<- pred10-pred00
    
  data1$d1<- pred11-pred01
  
  #Step 3: fit two regressions to model d0 and d1
  modeld0 <-gbm(d0~., data = data0[!is.na(data0$d0),-which(names(data0)%in%c("Time","Event","r0"))], 
                distribution="gaussian",n.trees = ntrees,interaction.depth=6, n.minobsinnode = round(959/100),bag.fraction = 0.8)
  
  modeld1 <-gbm(d1~., data = data1[!is.na(data1$d1),-which(names(data1)%in%c("Time","Event","r1"))], 
                distribution="gaussian",n.trees = ntrees,interaction.depth=6, n.minobsinnode = round(959/100),bag.fraction = 0.8)
  
  predictiond0<-predict(modeld0,testdat,n.trees=ntrees)
  predictiond1<-predict(modeld1,testdat,n.trees=ntrees)

  #step 4: fit propensity model of treatment assignment
  #random forest
  data$Treatment<-as.factor(data$Treatment)
  rf<-randomForest(Treatment~., data=data[,-which(names(data)%in%c("Time","Event"))],
                   ntree = ntrees, classwt = c(0.5, 0.5))
  pred<-predict(rf,testdat,type="prob")[,2]
  
  X.diff<-pred*predictiond0+(1-pred)*predictiond1

  return(list(diff=X.diff))
}



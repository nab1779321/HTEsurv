DNNsurv_HTE_T <- function(data=dat,testdat=testdat,time.interest=NULL){
  data<-data[complete.cases(data),]
  data1 <- data[which(data$Treatment==1),-c(which(names(data)=="Treatment"))]
  data0 <- data[which(data$Treatment==0),-c(which(names(data)=="Treatment"))]
  
  # prepare training data for treatment group =1
  x_dat1 = as.matrix(data1[, -which(colnames(data1) %in% c("Time","Event"))], nrow = nrow(data1))
  y_dat1 = as.matrix(data1[, c("Time","Event")], nrow = nrow(data1)) 
  
  # prepare training data for treatment group =0
  x_dat0 = as.matrix(data0[, -which(colnames(data0) %in% c("Time","Event"))], nrow = nrow(data0))
  y_dat0 = as.matrix(data0[, c("Time","Event")], nrow = nrow(data0))   
  
  # prepare prediction data
  x_val = as.matrix(testdat[, -which(colnames(testdat) %in% c("Time","Event","Treatment"))], nrow = nrow(testdat))  
  #y_val =as.matrix(data[, c("Time","Event")], nrow = nrow(data)) # notice: col1 is time, col2 is event
  
  # set up DNN parameters 
  num_nodes <- 30 # number of nodes per hidden layer
  string_activation <- "selu" # activation function
  num_l1 <- 0.1 # L1 penalty
  num_lr <- 0.01 # learning rate
  num_epoch <- 80 # number of epoches for optimization
  batch_size <- 50 # number of batch size for optimization
  
  #-------------------- Fit DNNsurv for treatment group=1 ----------------------#
  
  # Build the Keras model
  model1 <- keras_model_sequential() %>%
    layer_dense(units = num_nodes, activation = string_activation, input_shape = c(ncol(x_dat1)), kernel_regularizer = regularizer_l1(num_l1)) %>%
    layer_dense(units = num_nodes, activation = string_activation, kernel_regularizer = regularizer_l1(num_l1)) %>%
    layer_dense(units = 1, activation = string_activation)
  
  summary(model1)
  
  # Compile the Keras model
  model1 %>% compile(
    optimizer = optimizer_rmsprop(lr = num_lr),
    loss = loss_lik_efron,
    metrics = NULL)
  
  # Train the model 
  history1 <- model1 %>% keras::fit(x_dat1, y_dat1, epochs = num_epoch, batch_size = batch_size)
  plot(history1)
  
  #-------------------- Fit DNNsurv for treatment group=0 ----------------------#
  
  # Build the Keras model
  model0 <- keras_model_sequential() %>%
    layer_dense(units = num_nodes, activation = string_activation, input_shape = c(ncol(x_dat0)), kernel_regularizer = regularizer_l1(num_l1)) %>%
    layer_dense(units = num_nodes, activation = string_activation, kernel_regularizer = regularizer_l1(num_l1)) %>%
    layer_dense(units = 1, activation = string_activation)
  
  summary(model0)
  
  # Compile the Keras model
  model0 %>% compile(
    optimizer = optimizer_rmsprop(lr = num_lr),
    loss = loss_lik_efron,
    metrics = NULL)
  
  # Train the model 
  history0 <- model0 %>% keras::fit(x_dat0, y_dat0, epochs = num_epoch, batch_size = batch_size)
  plot(history0)
  
  #----------------- predict prognostic index in treatment group=1 --------------------#
  library(survival)
  y_val_pred1 <- model1 %>% predict(x_val)
  y_val_pred0 <- model0 %>% predict(x_val) 
  
  #---------------- predict survival probability --------------------#
  # baseline cumulative hazard for the model of treatment group=1
  y_dat_pred1 <- model1 %>% predict(x_dat1)
  baseline1 <- base_efron(y_test=y_dat1, y_test_pred=y_dat_pred1) # notice: col1 in y_test is time, col2 is event
  Lambda_t1 = sapply(time.interest, function(x){baseline1$cumhazard$hazard[sum(baseline1$cumhazard$time <= x)]  })
  
  # baseline cumulative hazard for the model of treatment group=0
  y_dat_pred0 <- model0 %>% predict(x_dat0)
  baseline0 <- base_efron(y_test=y_dat0, y_test_pred=y_dat_pred0) # notice: col1 in y_test is time, col2 is event
  Lambda_t0 = sapply(time.interest, function(x){baseline0$cumhazard$hazard[sum(baseline0$cumhazard$time <= x)]  })
  
  # predict survival probability
  S_x1 = exp(-1 * exp(y_val_pred1) %*% matrix(Lambda_t1, nrow = 1))
  S_x0 = exp(-1 * exp(y_val_pred0) %*% matrix(Lambda_t0, nrow = 1)) 
  
  #-------------- calculate ITE -----------------#
  pred0<-S_x0
  pred1<-S_x1
  DNNsurv.diff <- pred1-pred0
  return(list(diff=DNNsurv.diff,pred0=pred0,pred1=pred1))
  
}



DNNsurv_HTE_X <- function(data=dat,testdat=testdat,time.interest=NULL,ntrees=1000){
  data<-data[complete.cases(data),]
  data1 <- data[which(data$Treatment==1),-c(which(names(data)=="Treatment"))]
  data0 <- data[which(data$Treatment==0),-c(which(names(data)=="Treatment"))]
  
  # prepare training data for treatment group =1
  x_dat1 = as.matrix(data1[, -which(colnames(data1) %in% c("Time","Event"))], nrow = nrow(data1))
  y_dat1 = as.matrix(data1[, c("Time","Event")], nrow = nrow(data1)) 
  
  # prepare training data for treatment group =0
  x_dat0 = as.matrix(data0[, -which(colnames(data0) %in% c("Time","Event"))], nrow = nrow(data0))
  y_dat0 = as.matrix(data0[, c("Time","Event")], nrow = nrow(data0))   
  
  # prepare training data for all individuals
  x_dat_all = as.matrix(data[, -which(colnames(data) %in% c("Time","Event"))], nrow = nrow(data))
  y_dat_all = as.matrix(data[, c("Time","Event")], nrow = nrow(data))   
  
  # set up DNN parameters 
  num_nodes <- 30 # number of nodes per hidden layer
  string_activation <- "selu" # activation function
  num_l1 <- 0.1 # L1 penalty
  num_lr <- 0.01 # learning rate
  num_epoch <- 80 # number of epoches for optimization
  batch_size <- 50 # number of batch size for optimization
  
  # step 1: fit 2 DNNsurv models
  #-----Fit DNNsurv for treatment group=1 -----#
 
  # Build the Keras model
  model1 <- keras_model_sequential() %>%
    layer_dense(units = num_nodes, activation = string_activation, input_shape = c(ncol(x_dat1)), kernel_regularizer = regularizer_l1(num_l1)) %>%
    layer_dense(units = num_nodes, activation = string_activation, kernel_regularizer = regularizer_l1(num_l1)) %>%
    layer_dense(units = 1, activation = string_activation)
  
  #summary(model1)
  
  # Compile the Keras model
  model1 %>% compile(
    optimizer = optimizer_rmsprop(lr = num_lr),
    loss = loss_lik_efron,
    metrics = NULL)
  
  # Train the model 
  history1 <- model1 %>% keras::fit(x_dat1, y_dat1, epochs = num_epoch, batch_size = batch_size)
  
  #---- Fit DNNsurv for treatment group=0 ----#
  
  # Build the Keras model
  model0 <- keras_model_sequential() %>%
    layer_dense(units = num_nodes, activation = string_activation, input_shape = c(ncol(x_dat0)), kernel_regularizer = regularizer_l1(num_l1)) %>%
    layer_dense(units = num_nodes, activation = string_activation, kernel_regularizer = regularizer_l1(num_l1)) %>%
    layer_dense(units = 1, activation = string_activation)
  
  #summary(model0)
  
  # Compile the Keras model
  model0 %>% compile(
    optimizer = optimizer_rmsprop(lr = num_lr),
    loss = loss_lik_efron,
    metrics = NULL)
  
  # Train the model 
  history0 <- model0 %>% keras::fit(x_dat0, y_dat0, epochs = num_epoch, batch_size = batch_size)
  
  # Step 2: calculate D1 and D2
  # predict prognostic index
  library(survival)
  y_dat_pred01 <- model0 %>% predict(x_dat1) # trt data, plug in ctrl model
  y_dat_pred11 <- model1 %>% predict(x_dat1) # trt data, plug in trt model
  y_dat_pred00 <- model0 %>% predict(x_dat0) # ctrl data, plug in ctrl model
  y_dat_pred10 <- model1 %>% predict(x_dat0) # ctrl data, plug in trt model 
  
  # baseline cumulative hazard for the model of treatment group=1
  baseline1 <- base_efron(y_test=y_dat1, y_test_pred=y_dat_pred11) # notice: col1 in y_test is time, col2 is event
  Lambda_t1 = sapply(time.interest, function(x){baseline1$cumhazard$hazard[sum(baseline1$cumhazard$time <= x)]  })
  
  # baseline cumulative hazard for the model of treatment group=0
  baseline0 <- base_efron(y_test=y_dat0, y_test_pred=y_dat_pred00) # notice: col1 in y_test is time, col2 is event
  Lambda_t0 = sapply(time.interest, function(x){baseline0$cumhazard$hazard[sum(baseline0$cumhazard$time <= x)]  })
  
  # predict survival probability
  S_x01 = exp(-1 * exp(y_dat_pred01) %*% matrix(Lambda_t0, nrow = 1)) # trt data, put in ctrl group, use baseline cumulative hazard from ctrl group
  S_x11 = exp(-1 * exp(y_dat_pred11) %*% matrix(Lambda_t1, nrow = 1)) # trt data, put in trt group, use baseline cumulative hazard from trt group
  S_x00 = exp(-1 * exp(y_dat_pred00) %*% matrix(Lambda_t0, nrow = 1)) # ctrl data, put in ctrl group, use baseline cumulative hazard from ctrl group
  S_x10 = exp(-1 * exp(y_dat_pred10) %*% matrix(Lambda_t1, nrow = 1)) # ctrl data, put in trt group, use baseline cumulative hazard from trt group
  
  # D1=mu1(1)-mu0(1)
  data1$d1<- S_x11-S_x01
    
  # D0=mu1(0)-mu0(0)
  data0$d0<- S_x10-S_x00
  
  #Step 3: fit two regressions to model d0 and d1
  modeld0 <-gbm(d0~., data = data0[,-which(names(data0)%in%c("Time","Event","r0","log_r0"))], 
                distribution="gaussian",n.trees = ntrees,interaction.depth=6, n.minobsinnode = round(959/100),bag.fraction = 0.5)
  
  modeld1 <-gbm(d1~., data = data1[,-which(names(data1)%in%c("Time","Event","r1","log_r1"))], 
                distribution="gaussian",n.trees = ntrees,interaction.depth=6, n.minobsinnode = round(959/100),bag.fraction = 0.5)
  
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


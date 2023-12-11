import os
dirc='Your working directory'
os.chdir(dirc)
    
import pyreadr
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense
import tensorflow as tf
import tensorflow.keras.backend as tfkb
from tensorflow.keras import regularizers
from fun_util import loss_lik_efron, base_efron
import pandas as pd
import matplotlib.pylab as plt
from tensorflow.keras.callbacks import ModelCheckpoint
from sklearn import ensemble
import tensorflow as tf
import tensorflow.keras.backend as tfkb
import numpy as np
import multiprocessing
tf.config.run_functions_eagerly(True)

# specify the time of interest, train and test data, number of simulation
time_interest = np.array(15)
name_of_train_data='dat_all_weib_balanced_mediantt_n1000p10_test10000_scenario1.RData'
name_of_test_data='testDat10000_weib_mediantt_balanced_setting1_dataframe.RData'
num_sim = 100
output_name = "Xlearner_result_weib_balanced_mediantt_scenario1.csv"

TrainData = pyreadr.read_r(name_of_train_data)['dat_all']
TestData = pyreadr.read_r(name_of_test_data)['testdat']

TrainData.dropna(axis = 0, how = 'any', inplace = True)
TestData.dropna(axis = 0, how = 'any', inplace = True)

def Xlearner(sim_index):
    
    data = TrainData[(TrainData['sim']==sim_index)]
    data1 = TrainData[(TrainData['Treatment']==1) & (TrainData['sim']==sim_index)]
    data0 = TrainData[(TrainData['Treatment']==0) & (TrainData['sim']==sim_index)]
    
    # prepare training data for treatment group=1
    x_dat1 = tf.convert_to_tensor(data1.drop(['Treatment',"Time","Event",'sim'],axis=1))
    y_dat1 = tf.convert_to_tensor(data1[["Time","Event"]])
    
    # prepare training data for treatment group=0
    x_dat0 = tf.convert_to_tensor(data0.drop(['Treatment',"Time","Event",'sim'],axis=1))
    y_dat0 = tf.convert_to_tensor(data0[["Time","Event"]])
    
    # prepare prediction data
    x_val = tf.convert_to_tensor(TestData.drop(['Treatment',
                                                "Time",
                                                "Event",
                                                'true.diff',
                                                'true1',
                                                'true0'],axis=1))
    
    ### set up DNN parameters ###
    num_nodes = 30             # number of nodes per hidden layer
    #string_activation = "selu" # activation function
    num_l1 = 0.1               # L1 penalty
    num_lr = 0.01              # learning rate
    num_epoch = 80             # number of epoches for optimization
    batch_size = 50            # number of batch size for optimization
    
    #-------------------- Fit DNNsurv for treatment group=1 ----------------------#
    #clear the previous model
    tfkb.clear_session()
    
    # define the keras model
    model1 = Sequential()
    model1.add(Dense(num_nodes, input_dim=x_dat1.shape[1], activation='selu',kernel_regularizer=regularizers.l1(num_l1)))
    model1.add(Dense(num_nodes, activation='selu',kernel_regularizer=regularizers.l1(num_l1)))
    model1.add(Dense(1, activation='selu'))
    
    # compile the keras model
    model1.compile(loss=loss_lik_efron, optimizer='RMSprop', metrics=None)
    tfkb.set_value(model1.optimizer.learning_rate, num_lr)
    
    # save the history of training
    #filepath="model_{epoch:02d}-{loss:.2f}.hdf5"
    #checkpoint = ModelCheckpoint(os.path.join(dirc,filepath),
    #                             save_weights_only=True,
    #                             monitor='loss',
    #                             save_freq="epoch")
    
    # fit the keras model on the training set
    
    
    history1=model1.fit(x_dat1, y_dat1, 
                        epochs=num_epoch, 
                        batch_size=batch_size,
                        verbose = 0)
    
    #                   callbacks=[checkpoint])
    
    #    plt.plot(history1.history['loss'])
    #    plt.xlabel('epoch')
    #    plt.ylabel('loss')
    #    plt.show()
    #    
    #-------------------- Fit DNNsurv for treatment group=0 ----------------------#
    #clear the previous model
    tfkb.clear_session()
    
    # define the keras model
    model0 = Sequential()
    model0.add(Dense(num_nodes, input_dim=x_dat0.shape[1], activation='selu',kernel_regularizer=regularizers.l1(num_l1)))
    model0.add(Dense(num_nodes, activation='selu',kernel_regularizer=regularizers.l1(num_l1)))
    model0.add(Dense(1, activation='selu'))
    
    # compile the keras model
    model0.compile(loss=loss_lik_efron, optimizer='RMSprop', metrics=None)
    tfkb.set_value(model0.optimizer.learning_rate, num_lr)
    
    # save the history of training
    #filepath="model_{epoch:02d}-{loss:.2f}.hdf5"
    #checkpoint = ModelCheckpoint(os.path.join(dirc,filepath),
    #                             save_weights_only=True,
    #                             monitor='loss',
    #                             save_freq="epoch")
    
    # fit the keras model on the training set
    
    history0=model0.fit(x_dat0, y_dat0, 
                        epochs=num_epoch, 
                        batch_size=batch_size,
                        verbose = 0)
    #                   callbacks=[checkpoint])
    										
    tfkb.clear_session()
    
    ### step 2: calculate D1 and D2
    # predict prognostic index
    
    y_dat_pred01 = model0.predict(x_dat1) # trt data, plug in ctrl model
    y_dat_pred11 = model1.predict(x_dat1) # trt data, plug in trt model
    y_dat_pred00 = model0.predict(x_dat0) # ctrl data, plug in ctrl model
    y_dat_pred10 = model1.predict(x_dat0) # ctrl data, plug in trt model 
    
    # baseline cumulative hazard for the model of treatment group=1
    baseline1 = base_efron(y_test=y_dat1, y_test_pred=y_dat_pred11)
    Lambda_t1 = baseline1['cumhazard'][np.max(np.where(baseline1['time'] <= time_interest))]
    
    # baseline cumulative hazard for the model of treatment group=0
    baseline0 = base_efron(y_test=y_dat0, y_test_pred=y_dat_pred00)
    Lambda_t0 = baseline0['cumhazard'][np.max(np.where(baseline0['time'] <= time_interest))]
    
    # predict survival probability
    S_x01 = np.exp(-1 * np.outer(np.exp(y_dat_pred01),Lambda_t0))
    S_x11 = np.exp(-1 * np.outer(np.exp(y_dat_pred11),Lambda_t1))
    S_x00 = np.exp(-1 * np.outer(np.exp(y_dat_pred00),Lambda_t0))
    S_x10 = np.exp(-1 * np.outer(np.exp(y_dat_pred10),Lambda_t1))
    
    #np.save('pred1',S_x1)
    #np.save('pred0',S_x0)
    # D1=mu1(1)-mu0(1)
    d1 = S_x11-S_x01
    #r1 = S_x11/S_x01
    
    # D0=mu1(0)-mu0(0)
    d0 = S_x10-S_x00
    #r0 = S_x10/S_x00
    
    # Step 3: fit two regressions to model d0 and d1
    params = {'n_estimators': 1000,
              'max_depth': 6,
              'min_samples_leaf': round(959/100),
              'subsample': 0.5,
              'loss': 'ls'}
    
    reg0 = ensemble.GradientBoostingRegressor(**params)
    reg0.fit(data0.drop(['Treatment',"Time","Event",'sim'],axis=1), d0.reshape(-1))
    predictiond0 = reg0.predict(TestData.drop(['Treatment',
                                               "Time",
                                               "Event",
                                               'true.diff',
                                               'true1',
                                               'true0'],axis=1))
    
    reg1 = ensemble.GradientBoostingRegressor(**params)
    reg1.fit(data1.drop(['Treatment',"Time","Event",'sim'],axis=1), d1.reshape(-1))
    predictiond1 = reg1.predict(TestData.drop(['Treatment',
                                               "Time",
                                               "Event",
                                               'true.diff',
                                               'true1',
                                               'true0'],axis=1))    
    
    # step 4: fit propensity model of treatment assignment
    rf = ensemble.RandomForestClassifier(n_estimators=1000,
                                          min_samples_leaf=5,
                                          max_features=1/3)
    rf.fit(data.drop(['Treatment',"Time","Event",'sim'],axis=1),data['Treatment'])
    pred = rf.predict_proba(TestData.drop(['Treatment',
                                           "Time",
                                           "Event",
                                           'true.diff',
                                           'true1',
                                           'true0'],axis=1))
    
    preddiff = pred[:,0]*predictiond0 + pred[:,1]*predictiond1
    
    predtrue = np.concatenate([preddiff,
                               TestData['true.diff']]).reshape(-1,1)
    return predtrue

mymap = map(Xlearner,
            range(1,num_sim+1))
result = list(mymap)
result = np.squeeze(np.array(result)).T
np.savetxt(output_name, result, delimiter=",")
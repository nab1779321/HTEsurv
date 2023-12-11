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
output_name = 'Tlearner_result_weib_balanced_mediantt_scenario1.csv'

TrainData = pyreadr.read_r(name_of_train_data)['dat_all']
TestData = pyreadr.read_r(name_of_test_data)['testdat']

TrainData.dropna(axis = 0, how = 'any', inplace = True)
TestData.dropna(axis = 0, how = 'any', inplace = True)

def Tlearner(sim_index):

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
    
    y_val_pred1 = model1.predict(x_val)
    y_val_pred0 = model0.predict(x_val)
    										
    tfkb.clear_session()
    
    #---------------- predict survival probability --------------------#
    # baseline cumulative hazard for the model of treatment group=1
    y_dat_pred1 = model1.predict(x_dat1)
    baseline1 = base_efron(y_test=y_dat1, y_test_pred=y_dat_pred1)
    Lambda_t1 = baseline1['cumhazard'][np.max(np.where(baseline1['time'] <= time_interest))]
        
    # baseline cumulative hazard for the model of treatment group=0
    y_dat_pred0 = model0.predict(x_dat0)
    baseline0 = base_efron(y_test=y_dat0, y_test_pred=y_dat_pred0)
    Lambda_t0 = baseline0['cumhazard'][np.max(np.where(baseline0['time'] <= time_interest))]
    
    # predict survival probability
    S_x1 = np.exp(-1 * np.outer(np.exp(y_val_pred1),Lambda_t1))
    S_x0 = np.exp(-1 * np.outer(np.exp(y_val_pred0),Lambda_t0))
    
    #np.save('pred1',S_x1)
    #np.save('pred0',S_x0)
    
    predtrue = np.concatenate([(S_x1-S_x0).reshape(-1),
                               TestData['true.diff']]).reshape(-1,1)
    return predtrue

mymap = map(Tlearner,
            range(1,num_sim+1))
result = list(mymap)
result = np.squeeze(np.array(result)).T
np.savetxt(output_name, result, delimiter=",")
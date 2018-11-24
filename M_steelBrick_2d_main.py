# =================================
# Author: T.H.E. Gulikers
# Main file (Neural network model)
# Description:
# NN model of steel elementary brick, 2D, version 2
# =================================

#region import modules
import pandas as pd
import matplotlib.pyplot as plt
import os
from NN_helpers import *
os.environ["CUDA_VISIBLE_DEVICES"]="-1"   # don't run on GPU.. for some reason faster on my pc. comment to run on GPU
import tensorflow as tf
import time

tic = time.time()
#endregion

#region ======== SET PARAMETERS =========

# data
path        = 'C:/Users/tom/Documents/ANN_models/data_steelBrick_2D_1' # path were data is to be found
filename    = 'T_steelBrick'                                # if multiple files are to be imported, only provide common
                                                            # part of the name. Always have .txt as extension!
plotFileNames = ['T_steelBrick_L10_W10_D1_Sa450_Sb300_dS4',
                 'T_steelBrick_L10_W10_D1_Sa150_Sb450_dS4'] # list of file(s) to be plotted as result AND serve as
                                                            # test data

sigmaCap = 500.                                             # cut off data set above this value

# network properties
layer_def = [[10, 'relu'],[10, 'relu'], [None, 'tanh']]  # sequence of activation functions and nodes to use per layer
                                                          # final layer size is arbitrary: auto set to output layer size

# set optimizer type and parameters. To change these, check https://keras.io/optimizers/
# learn_rate, decay, beta1, eps = 0,0,0,0
optimFunc = keras.optimizers.Nadam(lr=0.008, schedule_decay=0.005, beta_1=0.95, beta_2=0.999, epsilon=1e-8)
lossFunc = 'mean_absolute_error'

# Sets weight initiation values of each layer. To change these, check https://keras.io/initializers/
kernel_initializer = 'lecun_normal'
bias_initializer = 'random_uniform'

# data parameters: set the maximum and minimum stress to which the normalisation of in and output data is performed
stressNorm_maxmin = np.array([[500., -500.]]).T
dstressNorm_maxmin = np.array([[15., -15.]]).T
strainNorm_maxmin = np.array([[.17, -.17]]).T
dstrainNorm_maxmin = np.array([[.005, -0.005]]).T

n_epochs, n_batch = 1000, 40000      # training parameters; epoch count and number of data points per batch
#endregion


#region import data
pathFileLst = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]    # total list of files in path
fileLst = []                                # filtered list of output files in path
for i in range(len(pathFileLst)):
    if ('.txt' in pathFileLst[i]) and (filename in pathFileLst[i]):
        fileLst.append(pathFileLst[i])

# build data arrays for training and plotting/testing
dataArr = pd.DataFrame()                                                    # extract linear part of curve
dataArr_plot, dataArr_plot_info, i = [], [], 0                              # initiate iterable list for plot data

plotFileLst = [f for f in fileLst if f.strip('.txt') in plotFileNames]      # determine files which are to be plotted

for file in fileLst:
    dataArrInstance = pd.read_csv(path+'/'+file, header=0)

    # get location where plasticity begins
    iCap = dataArrInstance.index.get_loc(((dataArrInstance['Sigma_1'] - sigmaCap) ** 2).idxmin())

    if file in plotFileLst:
        dataArr_plot_info.append(getSimProperties(file.strip('.txt')))
        dataArr_plot.append(dataArrInstance.iloc[:iCap])
    else:
        dataArr = dataArr.append(dataArrInstance.iloc[:iCap], ignore_index=True)

# titles / keys of the pandas data columns
keys_in     = dataArr.keys()[:6]            # keys that represent the inputs to the network
keys_inc    = dataArr.keys()[4:6]           # keys that represent the increment of stress or strain
keys_out    = dataArr.keys()[6:]            # keys that represent the outputs of the network
keys_total  = keys_in.append(keys_out)      # all keys

# filter required data columns from the total data set
for key in dataArr.keys():
    if key not in keys_total:
        dataArr.drop(columns=key)
        for i in range(len(dataArr_plot)):      # this is a list of pandas objects, in contrast to dataArr (single pandas object)
            dataArr_plot[i] = dataArr_plot[i].drop(columns=key)
#endregion


#region process data to normalized training and testing sets

# save max, min in dataframe to normalise and de-normalise
normparams = pd.DataFrame(np.concatenate((stressNorm_maxmin, stressNorm_maxmin, strainNorm_maxmin, strainNorm_maxmin,
                            dstressNorm_maxmin, dstressNorm_maxmin, dstrainNorm_maxmin, dstrainNorm_maxmin), axis=1),
                          index=['max', 'min'], columns=dataArr.keys())
dataArr_norm = normalise_data(dataArr, normparams)                             # normalise data
dataArr_plot_norm = pd.DataFrame()
for i in range(len(dataArr_plot)):
    dataArr_plot_norm = dataArr_plot_norm.append(normalise_data(dataArr_plot[i], normparams), ignore_index=True)  # normalise plot/test data
# dataArr_norm = dataArr_norm.sample(frac=1).reset_index(drop=True)                           # shuffle columns


# split into train and test data
train, test = dataArr_norm, dataArr_plot_norm

# creation of separate data arrays for in and outputs
train_in, train_out, test_in, test_out = train[keys_in], train[keys_out], test[keys_in], test[keys_out]
n_inputs, n_outputs = train_in.shape[1], train_out.shape[1]         # extract amount of in and outputs of the network
#endregion

#region define and train the network
toc = time.time()
print('Preprocessing finished at T = ', toc-tic)
# assemble network
network_geometry = init_network(layer_def, n_inputs, n_outputs, kernel_initializer, bias_initializer)  # define NN structure
network = compile_network(network_geometry, set_optim=optimFunc, set_floss=lossFunc)    # compile NN

# training process
history = network.fit(train_in, train_out, epochs=n_epochs, batch_size=n_batch, validation_data=(test_in, test_out))
toc = time.time()
print('NN training finished at T = ', toc-tic)
#endregion

#region evaluate network performance
score_train = np.array(network.evaluate(train_in, train_out, batch_size=n_batch))[1]
score_test = np.array( network.evaluate(test_in, test_out, batch_size=n_batch))[1]
print("Training accuracy: ", 100.*score_train, " %")
print("Testing accuracy: ", 100.*score_test, " %")

# summarize history for accuracy
sample_rate = 20
plt.plot(np.arange(1,n_epochs+1,sample_rate), history.history['acc'][::sample_rate], lw=1)
plt.plot(np.arange(1,n_epochs+1,sample_rate), history.history['val_acc'][::sample_rate], lw=1)
plt.title('model accuracy')
plt.ylim([0,1])
plt.ylabel('accuracy')
plt.xlabel('epoch')
plt.legend(['train', 'test']) #loc='upper left'
plt.show()

plot_network(network)
#endregion

quit()
for i in range(len(dataArr_plot)):

    #region reconstruct stress-strain based on NN

    # reconstructed curve initiation: dataArr_recon will contain neural network predicted values
    sigMin, sigMax, steps = [0., 0.], dataArr_plot_info[i][3:5], 80        # set reconstruction min and max stresses
    dataArr_recon = pd.DataFrame(0., index=[0], columns=dataArr.keys())    # initiate zeros reconstruction array

    # initiate and compute (non-normalised) increments and use them to construct the first entry/increment of the array
    sig_delta0 = np.zeros((1, len(sigMin)))
    for j in range(len(sigMin)):
        sig_delta0[0, j] = (sigMax[j]*1.1 - sigMin[j]) / steps
    dataArr_recon[keys_inc] = sig_delta0

    # normalise first entry of reconstruction array
    dataArr_recon[keys_inc] = normalise_data(dataArr_recon, normparams)[keys_inc]

    # loop to construct stress-strain curves with established neural network
    for j in range(len(np.arange(steps))):
        stateCurrent = pd.DataFrame(dataArr_recon.loc[j:j+1]).reset_index(drop=True)      # grab previous state
        stateUpdate = NNpredict(network, stateCurrent[keys_in], normparams, keys_in, keys_inc, keys_out)
        dataArr_recon = dataArr_recon.append(stateUpdate.join(stateCurrent[keys_inc.append(keys_out)]), ignore_index=True)

    dataArr_recon = normalise_data(dataArr_recon, normparams, deNorm = True)
    toc = time.time()

    #endregion

    #region plot stress-strain curves
    tableau = load_plotColors()

    sig1, sig2 = sigMax
    ticks_strain = np.linspace(0., 0.12,5)
    ticks_stress = np.linspace(0., 5e2, 5)

    plt.figure(figsize=(15,4.5))

    ax = plt.subplot(121)
    plt.title('X-direction')

    plt.plot(dataArr_plot[i]['Epsilon_1'], dataArr_plot[i]['Sigma_1'], label='Abaqus output', color=tableau[0], lw=2)
    plt.plot(dataArr_recon['Epsilon_1'], dataArr_recon['Sigma_1'], label='ANN output', color=tableau[5], lw=2)
    plt.grid(True, 'major', 'y', linestyle='--', alpha=.3)
    plt.xlim(-0.01, 0.125)
    plt.ylim(-50, 5e2)
    plotStyle(ax, ticks_strain, ticks_stress)


    ax = plt.subplot(122)
    plt.title('Y-direction')

    plt.plot(dataArr_plot[i]['Epsilon_2'], dataArr_plot[i]['Sigma_2'], label='Abaqus output', color=tableau[0], lw=2)
    plt.plot(dataArr_recon['Epsilon_2'], dataArr_recon['Sigma_2'], label='ANN output', color=tableau[5], lw=2)
    plt.grid(True, 'major', 'y', linestyle='--', alpha=.3)
    plt.xlim(-0.01, 0.125)
    plt.ylim(-30, 5e2)
    plotStyle(ax, ticks_strain, ticks_stress)

    plt.legend()

    print('Stress-strain curve {} reconstruction finished at T = '.format(i), toc-tic)
    plt.show()
    #endregion


# np.savetxt('steel'+str(n_epochs)+'.txt', dataArr_recon)
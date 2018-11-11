# =================================
# Author: T.H.E. Gulikers
# Main file (Neural network model)
# Description:
# NN model of general elementary brick, 2D, version 4
# Includes shear
# =================================

#region import modules
import os
from NN_helpers import *
os.environ["CUDA_VISIBLE_DEVICES"]="-1"   # don't run on GPU.. for some reason faster on my pc. comment to run on GPU
import time

tic = time.time()
#endregion

#region ======== SET PARAMETERS =========

# data
path        = 'C:/Users/tom/Documents/ANN_models/data_compBrick_2D_2/alldata_3'           # path were data is to be found
filename    = 'T_compBrick'                 # if multiple files are to be imported, only provide common
                                            # part of the name.

modifiers = []                              # modify data set. 'Damage' adds a damage parameter,
                                            # 'history1' adds the historic stress/strain state at step -1
                                            # 'history2' adds the historic stress/strain state at step -2

# Abaqus UMAT subroutine stuff
umat_fname = 'umat_test.for'                # file name to be given to the UMAT subroutine
material_name = 'compBrick_2_2_3'           # material name for UMAT subroutine

# list of file(s) to be plotted as result AND serve as test data
plotFileNames = ['T_compBrick_L10_W10_D1_Ea13333_Eb20000_Es6666_dE5',
                 'T_compBrick_L10_W10_D1_Ea6666_Eb13333_Es20000_dE15',
                 'T_compBrick_L10_W10_D1_Ea0_Eb0_Es20000_dE40']

EpsilonCap = .025                                             # [strain * 10^6] cut off data set above this value

# network properties
layer_def = [[60, 'relu'], [60, 'relu'], [None, 'tanh']]  # sequence of activation functions and nodes to use per layer
                                                                        # final layer size is arbitrary: auto set to output layer size

# set optimizer type and parameters. To change these, check https://keras.io/optimizers/
learn_rate, decay, eps = 8e-3, 5e-3, 1e-8
optimFunc = keras.optimizers.Nadam(lr=learn_rate, beta_1=0.95, beta_2=0.999, epsilon=eps, schedule_decay=decay)
lossFunc = 'mean_absolute_error'    # logcosh    # mean_absolute_error

# Sets weight initiation values of each layer. To change these, check https://keras.io/initializers/
kernel_initializer = 'lecun_normal' #'glorot_normal' #'random_uniform'
bias_initializer = 'random_uniform'

# data parameters: set the maximum and minimum stress to which the normalisation of in and output data is performed
stressNorm_maxmin = np.array([[1000., -1000.]]).T
dstressNorm_maxmin = np.array([[100., -100.]]).T
strainNorm_maxmin = np.array([[.025, -.025]]).T
dstrainNorm_maxmin = np.array([[.004, -0.004]]).T
damageParam_maxmin = np.array([[1., -1.]]).T

n_epochs, n_batch = 10000, 100000      # training parameters; epoch count and number of data points per batch
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
    dataArrInstance = processCols(pd.read_csv(path+'/'+file, header=0), modifiers)   # load data and add derived columns

    # get location where plasticity begins
    iCap = np.max(( dataArrInstance.index.get_loc(((dataArrInstance['Epsilon_1'] - EpsilonCap).round(3) ** 2).idxmin()),
                    dataArrInstance.index.get_loc(((dataArrInstance['Epsilon_2'] - EpsilonCap).round(3) ** 2).idxmin()),
                    dataArrInstance.index.get_loc(((dataArrInstance['Epsilon_12'] - EpsilonCap).round(3) ** 2).idxmin())))
    if iCap<=20: iCap=10000
    if file in plotFileLst:
        dataArr_plot_info.append(getSimProperties(file.strip('.txt'), strainControlled=True, dims=3))
        dataArr_plot_info[-1][3] = dataArr_plot_info[-1][3] / 1e6       # convert microstrain to strain
        dataArr_plot_info[-1][4] = dataArr_plot_info[-1][4] / 1e6
        dataArr_plot_info[-1][5] = dataArr_plot_info[-1][5] / 1e6

        dataArr_plot.append(dataArrInstance.iloc[:iCap])
    else:  # for leaving the testing data out of the training set
        dataArr = dataArr.append(dataArrInstance.iloc[:iCap], ignore_index=True)

# titles / keys of the pandas data columns
keys = dataArr.keys()
keys_out    = keys[[i for i in range(len(keys)) if 'dSigma' in keys[i]]]    # keys that represent the NN outputs
keys_inc    = keys[[i for i in range(len(keys)) if 'dEpsilon' in keys[i]]]  # keys that represent the NN increments
keys_in = keys[[i for i in range(len(keys)) if (keys[i] not in keys_out) and ('Sigma_' in keys[i] or 'Epsilon_' in keys[i])]]

# additional keys corresponding to history and damage parameter additions.
keys_in = keys_in.append(keys[[i for i in range(len(keys)) if '-1' in keys[i] and 'History1' in modifiers]])
keys_in = keys_in.append(keys[[i for i in range(len(keys)) if '-2' in keys[i] and 'History2' in modifiers]])
keys_in = keys_in.append(keys[[i for i in range(len(keys)) if 'Damage' in keys[i] and 'Damage' in modifiers]])

keys_total  = keys_in.append(keys_out)                                      # all keys

# filter required data columns from the total data set
for key in keys:
    if key not in keys_total:
        dataArr.drop(columns=key)
        for i in range(len(dataArr_plot)):      # this is a list of pandas objects, in contrast to dataArr (single pandas object)
            dataArr_plot[i] = dataArr_plot[i].drop(columns=key)
keys = dataArr.keys()   # renew key array after dropping columns
#endregion

#region process data to normalized training and testing sets

# save max, min in dataframe to normalise and de-normalise
normparams = pd.DataFrame(0., index=['max', 'min'], columns=keys)
for key in keys:
    if 'dSigma' in key: normparams[key] = dstressNorm_maxmin
    elif 'Sigma' in key:normparams[key] = stressNorm_maxmin
    elif 'dEpsilon' in key: normparams[key] = dstrainNorm_maxmin
    elif 'Epsilon' in key:normparams[key] = strainNorm_maxmin
    elif 'Damage' in key: normparams[key] = damageParam_maxmin
    elif 'Sigma-' in key: normparams[key] = stressNorm_maxmin
    elif 'Epsilon-' in key: normparams[key] = strainNorm_maxmin

dataArr_norm = normalise_data(dataArr, normparams)                             # normalise data
dataArr_plot_norm = pd.DataFrame()
for i in range(len(dataArr_plot)):
    dataArr_plot_norm = dataArr_plot_norm.append(normalise_data(dataArr_plot[i], normparams), ignore_index=True)  # normalise plot/test data

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

#region evaluate network by plotting the test data set
for i in range(len(dataArr_plot)):

    # region reconstruct stress-strain based on NN

    # reconstructed curve initiation: dataArr_recon will contain neural network predicted values
    epsMin, epsMax, steps = [0., 0., 0.], dataArr_plot_info[i][3:6], 100         # set reconstruction min and max strains
    dataArr_recon = pd.DataFrame(0., index=[0], columns=dataArr.keys())     # initiate zeros reconstruction array

    # initiate and compute (non-normalised) increments and use them to construct the first entry/increment of the array
    eps_delta0 = np.zeros((1, len(epsMin)))
    for j in range(len(epsMin)):
        eps_delta0[0, j] = (epsMax[j] - epsMin[j]) / steps
    dataArr_recon[keys_inc] = eps_delta0

    # first damage parameter values
    try: dataArr_recon.loc[0, ['DamageParameter_1', 'DamageParameter_2']] = 1.
    except: pass

    # normalise first entry of reconstruction array
    dataArr_recon[keys_inc] = normalise_data(dataArr_recon, normparams)[keys_inc]

    # loop to construct stress-strain curves with established neural network
    for j in range(len(np.arange(steps))):
        stateCurrent = pd.DataFrame(dataArr_recon.loc[j:j+1]).reset_index(drop=True)      # grab previous state
        stateUpdate = NNpredict(network, stateCurrent[keys_in], normparams, keys_in, keys_inc, keys_out)
        stateNew = pd.DataFrame(0., columns=stateCurrent.keys(), index=[0])        # empty array
        stateNew[stateUpdate.keys()] = stateUpdate
        stateNew[keys_inc.append(keys_out)] = stateCurrent[keys_inc.append(keys_out)]
        dataArr_recon = dataArr_recon.append(stateNew, ignore_index=True, sort=False)

    dataArr_recon = normalise_data(dataArr_recon, normparams, deNorm = True)
    toc = time.time()

    #endregion

    #region plot stress-strain curves
    tableau = load_plotColors()
    eps1, eps2, eps12 = epsMax
    ticks_strain = np.linspace(0., 0.02,5)
    ticks_stress = np.linspace(0., 7e2, 5)
    
    plt.figure(figsize=(15,5))

    ax = plt.subplot(131)
    plt.title('X-direction')
    
    plt.plot(dataArr_plot[i]['Epsilon_1'], dataArr_plot[i]['Sigma_1'], label='Abaqus output', color=tableau[0], lw=2)
    plt.plot(dataArr_recon['Epsilon_1'], dataArr_recon['Sigma_1'], label='NN output', color=tableau[5], lw=2)
    plt.grid(True, 'major', 'y', linestyle='--', alpha=.3)
    plt.xlim(-0.001, 0.022)
    plt.ylim(-50, 7e2)
    plotStyle(ax, ticks_strain, ticks_stress)

    ax = plt.subplot(132)
    plt.title('Y-direction')

    plt.plot(dataArr_plot[i]['Epsilon_2'], dataArr_plot[i]['Sigma_2'], label='Abaqus output', color=tableau[0], lw=2)
    plt.plot(dataArr_recon['Epsilon_2'], dataArr_recon['Sigma_2'], label='NN output', color=tableau[5], lw=2)
    plt.grid(True, 'major', 'y', linestyle='--', alpha=.3)
    plt.xlim(-0.001, 0.022)
    plt.ylim(-50, 7e2)
    plotStyle(ax, ticks_strain, ticks_stress)

    ax = plt.subplot(133)
    plt.title('XY-direction')

    plt.plot(dataArr_plot[i]['Epsilon_12'], dataArr_plot[i]['Sigma_12'], label='Abaqus output', color=tableau[0], lw=2)
    plt.plot(dataArr_recon['Epsilon_12'], dataArr_recon['Sigma_12'], label='NN output', color=tableau[5], lw=2)
    plt.grid(True, 'major', 'y', linestyle='--', alpha=.3)
    plt.xlim(-0.001, 0.022)
    plt.ylim(-50, 7e2)
    plotStyle(ax, ticks_strain, ticks_stress)

    plt.suptitle('stress-strain curves for loading ratio {} / {} / {} '.format(eps1, eps2, eps12))
    print('Stress-strain curve {} reconstruction finished at T = '.format(i+1), toc-tic)
    plt.legend()
    plt.show()
    #endregion
#endregion
print(learn_rate, decay, eps)
print(layer_def)


# if 'y' in input('make UMAT? (y/n)'):
#     createUMAT(umat_fname, material_name, network, normparams, 3)
#     toc = time.time()
#     print('UMAT constructed with file name: "{}" at T = '.format(umat_fname,i + 1), toc - tic)
#
# for i in range(len(dataArr_plot)):
#
#     tableau = load_plotColors()
#     eps1, eps2, eps12 = epsMax
#     ticks_strain = np.linspace(0., 0.02, 5)
#     ticks_stress = np.linspace(0., 7e2, 5)
#
#     plt.figure(figsize=(15, 5))
#
#     ax = plt.subplot(131)
#     plt.title('X-direction')
#
#     plt.plot(dataArr_plot[i]['Epsilon_1'], dataArr_plot[i]['Sigma_1'], label='Abaqus output', color=tableau[0], lw=2)
#     plt.grid(True, 'major', 'y', linestyle='--', alpha=.3)
#     plt.xlim(-0.001, 0.022)
#     plt.ylim(-50, 7e2)
#     plotStyle(ax, ticks_strain, ticks_stress)
#
#     ax = plt.subplot(132)
#     plt.title('Y-direction')
#
#     plt.plot(dataArr_plot[i]['Epsilon_2'], dataArr_plot[i]['Sigma_2'], label='Abaqus output', color=tableau[0], lw=2)
#     plt.grid(True, 'major', 'y', linestyle='--', alpha=.3)
#     plt.xlim(-0.001, 0.022)
#     plt.ylim(-50, 7e2)
#     plotStyle(ax, ticks_strain, ticks_stress)
#
#     ax = plt.subplot(133)
#     plt.title('XY-direction')
#
#     plt.plot(dataArr_plot[i]['Epsilon_12'], dataArr_plot[i]['Sigma_12'], label='Abaqus output', color=tableau[0], lw=2)
#     plt.grid(True, 'major', 'y', linestyle='--', alpha=.3)
#     plt.xlim(-0.001, 0.022)
#     plt.ylim(-50, 7e2)
#     plotStyle(ax, ticks_strain, ticks_stress)
#
#     plt.suptitle('stress-strain curves for loading ratio {} / {} / {} '.format(eps1, eps2, eps12))
#     print('Stress-strain curve {} reconstruction finished at T = '.format(i + 1), toc - tic)
#     plt.legend()
#     plt.show()
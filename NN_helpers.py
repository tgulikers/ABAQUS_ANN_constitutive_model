
#region import modules
import numpy as np
import matplotlib.pyplot as plt
from keras.models import Model, Sequential
import keras
from keras.layers import Dense, Input, Concatenate, Lambda
from keras.utils import plot_model
import tensorflow as tf
import os
from tensorflow import set_random_seed
import pandas as pd
import random as rn
import warnings

np.set_printoptions(threshold=np.nan)
#endregion


# =---------------------------------------------------------------------------------------------------
# Function: processCols
# Description: adds extra columns to the dataset, based on existing columns
# inputs:
# dataArr - pandas dataframe
# modifiers    - 'history1' copies Sigma and Epsilon to a new column, shifting it by 1 index down
#              - 'history2' copies Sigma and Epsilon to a new column, shifting it by 2 indices down
#              - 'damage' adds a quantification of the damage state at each increment
# outputs:
# dataframe with additional columns
# =---------------------------------------------------------------------------------------------------
def processCols(dataArr, modifiers = ''):

    if 'History1' in modifiers:
        for i in [int(key.strip('dEpsilon_')) for key in dataArr.keys() if 'dEpsilon' in key]: # find dimension (1D, 2D or 3D)
            i = str(i)  # enumerate starting from 1 and convert to string
            dataArr.loc[1:, 'Sigma-1_'+i] = dataArr.loc[:(len(dataArr['Sigma_1'])-2), 'Sigma_' + i].values
            dataArr.loc[1:, 'Epsilon-1_'+i] = dataArr.loc[:(len(dataArr['Epsilon_1'])-2), 'Epsilon_'+i].values
        dataArr.fillna(value=0., inplace=True)

    if 'History2' in modifiers:
        for i in range(1, len([i for i in dataArr.keys() if 'dEpsilon' in i])+1):  # find dimension (1D, 2D or 3D)
            i = str(i)  # enumerate starting from 1 and convert to string
            dataArr.loc[2:, 'Sigma-2_' + i] = dataArr.loc[:(len(dataArr['Sigma_1']) - 3), 'Sigma_' + i].values
            dataArr.loc[2:, 'Epsilon-2_' + i] = dataArr.loc[:(len(dataArr['Epsilon_1']) - 3), 'Epsilon_' + i].values
        dataArr.fillna(value=0., inplace=True)

    if 'Damage' in modifiers:     # ratio of current stress/strain slope compared to initial one
        for i in range(1, len([i for i in dataArr.keys() if 'dEpsilon' in i])+1):        # find dimension (1D, 2D or 3D)
            i = str(i)        # enumerate starting from 1 and convert to string
            dataArr['DamageParameter_'+i] = dataArr['Sigma_'+i] / dataArr['Epsilon_'+i] * \
                                           dataArr['Epsilon_'+i][1] / dataArr['Sigma_'+i][1]
        dataArr.fillna(value=1., inplace=True)
    return dataArr


# =---------------------------------------------------------------------------------------------------
# Function: normalise_data
# Description: normalises data between -1 and 1 for use with a tanh output function
# inputs:
# dataArr: pandas dataframe or numpy data array
# normalisationParams: pandas dataframe with amount of columns equal to dataArr and 2 rows labeled max, min
# outputs:
# dataframe with normalised columns
# =---------------------------------------------------------------------------------------------------
def normalise_data(dataArr, normalisationParams, deNorm=False):
    if not deNorm:      # loop for regular normalisation
        return 2.*(dataArr-normalisationParams.loc['min'])/ \
                (normalisationParams.loc['max'] - normalisationParams.loc['min'])-1.
    else:               # loop for denormalisation
        return (dataArr+1.)/2. * (normalisationParams.loc['max'] - normalisationParams.loc['min']) + \
                            normalisationParams.loc['min']


# =---------------------------------------------------------------------------------------------------
# Function: init_network
# Description: defines a neural network
# inputs:
# layer_def. List object in format [[#neurons 1, 'activation 1'], [[#neurons 2, 'activation 2']], ...]
# n_in, n_out: integer equal to number of network inputs, outputs
# act_out: activation of output layer. Standard is sigmoid
# outputs:
# model
# =---------------------------------------------------------------------------------------------------
def init_network(layer_def, n_in, n_out, initKernel, initBias, regularizer = None):
    model = Sequential()
    for n_layer in range(len(layer_def)-1):
        # first layer requires 'input_dim' argument; the if-statement enforces this.
        if n_layer==0:
            model.add(Dense(layer_def[n_layer][0], input_dim=n_in, kernel_initializer=initKernel,
                                       bias_initializer=initBias, activation=layer_def[n_layer][1],
                                       name='HiddenLayer_{}'.format(n_layer+1), bias_regularizer=regularizer))
            # model.add(keras.layers.Dropout(0.1))
        else:
            model.add(Dense(layer_def[n_layer][0], activation=layer_def[n_layer][1], kernel_initializer=initKernel,
                              bias_initializer=initBias, name='HiddenLayer_{}'.format(n_layer+1), bias_regularizer=regularizer))
        n = layer_def[n_layer][0]
    model.add(Dense(n_out, input_dim=n, activation=layer_def[n_layer+1][1], name='OutputLayer',
                    kernel_initializer=initKernel, bias_initializer=initBias))
    return model


# =---------------------------------------------------------------------------------------------------
# Function: compile_network
# Description: compiles a neural network model object to an object capable of being optimised
# inputs:
# NN model description
# set_optim: provide the optimiser to use in the training process
# set_floss: description of the loss function to be used to compute optimal steps
# set_metrics: leave as it is
# outputs:
# full NN model
# =---------------------------------------------------------------------------------------------------
def compile_network(model, set_optim = 'adam', set_floss='mean_square_error', set_metrics = ['accuracy']):
    model.compile(loss=set_floss, optimizer=set_optim, metrics=set_metrics)
    return model


# =---------------------------------------------------------------------------------------------------
# Function: NNpredict
# Description: performs NN computation of the trained network
# inputs:
# NN model description
# dataArr_in: Data set containing the required inputs to the network
# normalisationParams: pandas dataframe with amount of columns equal to dataArr and 2 rows labeled max, min
# keys_in, keys_out, keys_inc
# computeNextState - If false, NNpredict does not post-process the NN results
# outputs:
# pandas dataframe with (non-normalised!) updated stress-strain state
# =---------------------------------------------------------------------------------------------------
def NNpredict(model, dataArr_in, normalisationParams, keys_in, keys_inc, keys_out):
    if type(dataArr_in)!= pd.DataFrame: raise TypeError("Please provide a pandas dataframe as input!")

    # output in terms of denormalised strain increments
    delta_out = normalise_data(pd.DataFrame(model.predict(dataArr_in), columns=keys_out), normalisationParams[keys_out],
                               deNorm=True).reset_index(drop=True)

    # get input keys excl. increments
    keys_change = keys_in[[not bol for bol in -1 * np.isin(keys_in, keys_inc)]]

    # Initiate dataframe for final output
    dataArr_out = pd.DataFrame(0., index=[0], columns=normalisationParams.keys())
    dataArr_out[keys_in] = normalise_data(dataArr_in[keys_in], normalisationParams[keys_in], deNorm = True)
    dataArr_out[keys_out] = delta_out


    # construct dataframe of delta sigmas and delta epsilons
    delta = pd.DataFrame(0., index=[0], columns=normalisationParams.keys())
    delta[keys_inc] = dataArr_out[keys_inc]
    delta[keys_out] = delta_out

    #region update damage parameter
    # damageParameter_new = damageParameter_old * (epsilon_i / sigma_i) / (epsilon_i+1 / sigma_i+1)
    damageKeyLst = [key for key in dataArr_out.keys() if 'Damage' in key]
    for i_key in range(len(damageKeyLst)):
        ratio_in    = dataArr_out['Epsilon_' + str(i_key+1)] / dataArr_out['Sigma_' + str(i_key+1)]
        ratio_out   = (dataArr_out['Epsilon_' + str(i_key+1)] + delta['dEpsilon_' + str(i_key+1)]) / \
                      (dataArr_out['Sigma_' + str(i_key+1)] + delta['dSigma_' + str(i_key+1)])

        if (ratio_in / ratio_out).values[0] >= 1.: delta.loc[:, damageKeyLst[i_key]] = 1.
        elif (ratio_in / ratio_out).values[0] <= 0.: delta.loc[:,damageKeyLst[i_key]] = 0.
        else: delta.loc[:, damageKeyLst[i_key]] = (ratio_in / ratio_out).values[0]

        delta.fillna(value=1., inplace=True)        # change NaN's to 0

    warnings.filterwarnings("ignore")                               # suppress annoying and invalid warning in next line
    if len(damageKeyLst)>0: dataArr_out.loc[:, damageKeyLst] = 0.   # reset damage variable to avoid summing them in the next step
    warnings.filterwarnings("default")                              # display warnings again
    #endregion

    # update stress/strain states in output dataframe
    for key in keys_change:
        if '-2' in key: delta[key] = dataArr_out[key.replace('-2','-1')]          # historic stress/strain state
        elif '-1' in key: delta[key] = dataArr_out[key.replace('-1','')]            # historic stress/strain state
        elif 'Sigma_' in key or 'Epsilon_' in key:
            delta[key] = dataArr_out[key] + dataArr_out['d'+key]   # current stress/strain state

    # update dataArr_out to new stress/ strain state
    dataArr_out[keys_change] = delta[keys_change]

    return normalise_data(dataArr_out[keys_change], normalisationParams[keys_change])


# =---------------------------------------------------------------------------------------------------
# Function: getOdbProperties
# Description: extracts dimensions and loading from file name.
# inputs:
# jobname - format 'name'_L##_W##_D##_Sa##_Sb##.txt
# strainControlled - set to true if the data is strain controlled
# dims - amount of stress/strain components. Set to 2 for biaxial, set to 3 if including shear
# outputs:
# list of floats with [length, width, depth, sigma1, sigma2, (sigma3)] (stress controlled) or
# [length, width, depth, epsilon1, epsilon2, (epsilon3)] (strain controlled)
# =---------------------------------------------------------------------------------------------------
def getSimProperties(jobname, strainControlled=False, dims = 2):
    jobname += '_' # add underscore at end such that the program recognises the end of sigma2

    # find key letters indicating length, width, depth, load1, load2, dload
    locL, locW, locD = jobname.find('L'), jobname.find('W'), jobname.find('D')
    if not strainControlled: load1, load2, dload = jobname.find('Sa')+1, jobname.find('Sb')+1, jobname.find('dS')+1
    else: load1, load2, dload = jobname.find('Ea') + 1, jobname.find('Eb') + 1, jobname.find('dE') + 1

    if dims==3:
        if strainControlled: load3 = jobname.find('Es')+1
        else: load3 = jobname.find('Ss')+1

        lst, j = ['', '', '', '', '', ''], 0
        for loc in [locL, locW, locD, load1, load2, load3]:
            current, i, num = jobname[loc + 1], 2, []     # initialise loop parameters. current (currently considered string entry; an integer),
                                                          # i (index of current), num (list of valid numbers)
            while current != '_':                         # find how many numbers there are after each indicator letter
                num.append(current)                       # save integer to temporary list 'num'
                current = jobname[loc + i]                # extract next entry (to be checked for type == int)
                i += 1
            for k in num:                                 # converts num with format ['1', '2', '3'] to [123]
                lst[j] = str(lst[j]) + k
            j += 1
    else:
        lst, j = ['', '', '', '', ''], 0
        for loc in [locL, locW, locD, load1, load2]:
            current, i, num = jobname[loc + 1], 2, []     # initialise loop parameters. current (currently considered string entry),
                                                          # i (index of current), num (list of valid numbers)
            while current != '_':                         # find how many numbers there are after each letter
                num.append(current)                       # save integer to temporary list 'num'
                current = jobname[loc + i]                # extract next entry (to be checked for type == int)
                i += 1
            for k in num:                                 # converts num with format ['1', '2', '3'] to [123]
                lst[j] = str(lst[j]) + k
            j += 1

    return [float(i) for i in lst]


# =---------------------------------------------------------------------------------------------------
# Function: getNetworkWeights
# Description: extracts weights of each node from a NN model.
# inputs:
# model - a Keras/Tensorflow neural network object
# outputs:
# W - 2-D array with shape: [weights_layer_1, weights_layer_2, ... , weights_layer_N]
# =---------------------------------------------------------------------------------------------------
def getNetworkWeights(model):
    W = np.array(np.zeros((len(model.layers))), dtype=object)
    for i_layer in range(len(model.layers)):
        W[i_layer] = model.layers[i_layer].get_weights()
    return W


# =---------------------------------------------------------------------------------------------------
# Function: stressStrainToForceDisplacement
# Description: converts stress and strain increments to forces and displacements for a 2D 4-node quadrilateral element
# inputs:
# dsigma - array of stress increments
# depsilon - array of strain increments
# dims - list with dimensions of brick [length x width x thickness]
# Non-normalised data!
# outputs:
# df - array of force increments corresponding to DoF's
# du - array of displacement increments corresponding to DoF's
# DoF convention: up and right is positive. Numbering nodes: 1 to 4 starting top left, then clockwise
# Each node has 2 DoF's, so the element has 8. node 1, x dir is DoF[0]. node 4, y dir = DoF[7]
# =---------------------------------------------------------------------------------------------------
def stressStrainToForceDisplacement(dsigma, depsilon, dims):
    # convert
    Fx, Fy, Sx, Sy = dsigma[0]*dims[1]*dims[2], dsigma[1]*dims[0]*dims[2], 0.*dims[1]*dims[2], 0.*dims[0]*dims[2]
    Ux, Uy, Gx, Gy = depsilon[0]*dims[0], depsilon[1]*dims[1], 0.*dims[1], 0.*dims[0]   # G is shear deformation

    # create force vectors and displacement vectors
    # df := [f1x, f1y, f2x, f2y, f3x, f3y, f4x, f4y]
    df = np.array([Sx-Fx, -Sy+Fy, Sx+Fx, Sy+Fy, -Sx+Fx, Sy-Fy, -Sx-Fx, -Sy-Fy ])/2.

    # du := [f1x, f1y, f2x, f2y, f3x, f3y, f4x, f4y]
    du = np.array([Gx-Ux, -Gy+Uy, Gx+Ux, Gy+Uy, -Gx+Ux, Gy-Uy, -Gx-Ux, -Gy-Uy ])
    return df, du


# =---------------------------------------------------------------------------------------------------
# Function: computeKfromNN
# Description: Computes the stiffness matrix at a given stress-strain state, based on a NN model
# inputs:
# model - a Keras/Tensorflow neural network object
# state - data array with 1 set of inputs, defining the state and increment
# Normalised data!
# outputs:
# K - the local material stiffness matrix
# =---------------------------------------------------------------------------------------------------
def computeKfromNN(model, state, normparams, keys_in, network_geometry, n_dim):

    # find all required ingredients
    S_sigma = normparams.loc['max', 'Sigma_1']     # normalize ratio to get -1 < sigma_NN < 1 by: sigma_NN = sigma/S_sigma
    S_epsilon = normparams.loc['max', 'Epsilon_1'] # See above: -1 < (epsilon_NN = epsilon/S_epsilon) < 1
    beta = 1.                                      # tanh parameter. not used
    NB = model.layers[1].input_shape[1]            # length of hidden layer 1
    NC = model.layers[2].input_shape[1]            # length of hidden layer 2
    sigma_NN_i1 = model.predict(state[keys_in])[0] # network output (normalised)

    # weight matrices
    W = getNetworkWeights(model)
    w_BA = np.matrix(W[0][0]).T
    w_CB = np.matrix(W[1][0]).T
    w_DC = np.matrix(W[2][0]).T

    # node values for hidden layers 1 and 2
    B = networkOutputPartial(model, state, keys_in, 1)
    C = networkOutputPartial(model, state, keys_in, 2)

    # calculation of stiffness matrix according to HASHASH (2004) for a [relu, relu, tanh] network
    if network_geometry == ['relu', 'relu', 'tanh']:
        K = np.zeros((n_dim,n_dim))
        for i in range(n_dim):
            for j in range(n_dim):
                Part0 = S_sigma / S_epsilon * beta
                PartAB = 0.
                for k in range(NC):
                    PartA = ((1.-(sigma_NN_i1[i])**2)*w_DC[i, k])
                    PartB = 0.
                    for l in range(NB):
                        PartB += (min(np.ceil(max(C[k],0.)),1.) *w_CB[k, l])* (min(np.ceil(max(B[l],0.)),1.)*w_BA[l, j])
                    PartAB += PartA * PartB
                K[i,j] = Part0*PartAB
    # calculation of stiffness matrix for a [tanh, tanh, tanh] network
    elif network_geometry == ['tanh', 'tanh', 'tanh']:
        K = np.zeros((n_dim,n_dim))
        for i in range(n_dim):
            for j in range(n_dim):
                Part0 = S_sigma / S_epsilon * beta
                PartAB = 0.
                for k in range(NC):
                    PartA = ((1.-(sigma_NN_i1[i])**2)*w_DC[i, k])
                    PartB = 0.
                    for l in range(NB):
                        PartB += ((1.-(C[k])**2)*w_CB[k,l])* ((1.-(B[l])**2)*w_BA[l,j])
                    PartAB += PartA * PartB
                K[i,j] = Part0*PartAB
    # calculation of stiffness matrix for a [relu, relu, relu, tanh] network. DOUBLE-CHECK DERIVATION! NOT VERIFIED!
    elif network_geometry == ['relu', 'relu', 'relu', 'tanh']:
        w_ED, D, ND = np.matrix(W[3][0]).T, networkOutputPartial(model, state, keys_in, 3), model.layers[3].input_shape[1]
        K = np.zeros((n_dim,n_dim))
        for i in range(n_dim):
            for j in range(n_dim):
                Part0 = S_sigma / S_epsilon * beta
                PartABC = 0.
                for k in range(ND):
                    PartA = ((1.-(sigma_NN_i1[i])**2)*w_ED[k, i])
                    PartB = 0.    # loop over sum of l = 1 to NB
                    PartBC = 0.
                    for l in range(NC):
                        PartB += (min(np.ceil(max(D[k],0.)),1.) *w_DC[k, l])
                        PartC = 0.
                        for m in range(NB):
                            PartC = (min(np.ceil(max(C[l],0.)),1.)*w_CB[l, m])* \
                                    (min(np.ceil(max(B[m],0.)),1.)*w_BA[m, j])
                        PartBC += PartB*PartC
                    PartABC += PartA * PartBC
                K[i,j] = Part0*PartABC
    else:
        raise ModuleNotFoundError('Network architecture not supported by stiffness matrix calculation module')
    return K


# =---------------------------------------------------------------------------------------------------
# Function: computeKfromNN_2
# Description: Computes the stiffness matrix at a given stress-strain state, based on a NN model
# inputs:
# model - a Keras/Tensorflow neural network object
# state - data array with 1 set of inputs, defining the state and increment
# Non-normalised data!
# outputs:
# K - the local material stiffness matrix
# =---------------------------------------------------------------------------------------------------
def computeKfromNN_2(model, state, normparams, keys_in, keys_inc, keys_out, n_dim):
    K = np.zeros((n_dim,n_dim))

    dStrain = state[keys_inc].values[0]
    strainStates = ((dStrain[0]*4, dStrain[1]/2.), (dStrain[0]/2., dStrain[1]*4))

    # compute network output to inputs given in strainstates
    b = np.zeros((n_dim, n_dim))
    A = np.zeros((n_dim, n_dim))
    for i_state in range(len(strainStates)):
        warnings.filterwarnings("ignore")                               # suppress invalid warning in next line
        state.loc[:, keys_inc] = strainStates[i_state]
        warnings.filterwarnings("default")                              # display warnings again
        A[i_state, :] = normalise_data(state[keys_inc], normparams[keys_inc], deNorm=True)
        b[i_state, :] = normalise_data(model.predict(state[keys_in])[0], normparams[keys_out], deNorm=True)

    # solve systems of equations
    for i in range(len(strainStates)):
        K[i,:] = np.linalg.solve(A, b[:, i])

    return K


# =------------------------------------------------------------------------------------------------------
# Function: networkOutputPartial
# Description: computes the values of nodes at a requested hidden layer <i_layer_out>
# inputs:
# model - a Keras/Tensorflow neural network object
# state - data array with 1 set of inputs, defining the state and increment
# layer_def - specification of the network activation functions
# normalised data!
# outputs:
#
# =------------------------------------------------------------------------------------------------------
def networkOutputPartial(model, state, keys_in, i_layer_out):
    Wb = getNetworkWeights(model)
    x = state[keys_in].values.reshape((len(keys_in), 1))
    for i_layer in range(i_layer_out):
        W, b = Wb[i_layer]
        W, b = np.matrix(W.T), b.reshape((len(b), 1))

        activation = W * x + b
        if 'relu' in str(model.layers[i_layer].activation): output = np.multiply(activation,(activation>0))
        elif 'tanh' in str(model.layers[i_layer].activation): output = np.tanh(activation)
        else: raise LookupError("No activation function found in layer ", i_layer)

        x = output
    return x


# =------------------------------------------------------------------------------------------------------
# Function: createUMAT
# Description: writes a UMAT user-subroutine in fortran 95 for abaqus to represent a neural network model
# only works for a network with 1, 2, 3 or 4 hidden layers (easily adjustable for more)
# network structure currently implemented: All layers relu, except final layer which is tanh
# inputs:
# fname - file name
# model - a Keras/Tensorflow neural network object
# normalisationParams - pandas dataframe with amount of columns equal to dataArr and 2 rows labeled max, min
# n_dim - the dimensionality of the problem. For example, if sigma_1, sigma_2 and sigma_12, then n_dim = 3
# solver - either 'secant' (default), linear or 'hashash'. The latter option is based on the work of hashash (2004)
# outputs:
# D - the local material stiffness matrix
# =------------------------------------------------------------------------------------------------------
def createUMAT(fname, matname, model, normalisationParams, n_dim, solver = 'secant'):

    W = getNetworkWeights(model)
    w_BA, b_BA = W[0][0].T, W[0][1].reshape(-1,1)
    w_final, b_final = W[-1][0].T, W[-1][1].reshape(-1,1)

    f = open(fname, 'w')

    #region frontmatter

    # include external modules
    f.write('!****************************************!\n'
            '!   Abaqus UMAT neural network interface !\n'
            '!   specific for 2D shell element        !\n'
            '!   created by Tom Gulikers, TU Delft    !\n'
            '!   Date: 01-Okt-2018                    !\n'
            '!                                        !\n'
            '!****************************************!\n')
    f.write('! force free-form Fortran \n !DIR$ FREEFORM\n\n')
    f.write('!------ include external modules --------------------------\n'+
            " include 'globals/parameter_module.f90' \n"+
            '!------------------------------------------------------\n\n')

    # initiate umat subroutine by standard header
    f.write('SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,&\n'+
            ' &DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP,&\n'+
            ' &PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,&\n'+
            ' &COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,&\n'+
            ' &KSPT, KSTEP, KINC)\n\n'
            '  ! load FNM modules\n  use parameter_module,       only: NDIM, DP, ZERO, ONE, SMALLNUM\n\n'+
            "  include 'aba_param.inc'\n\n"+
            '  CHARACTER(len=8) :: CMNAME\n\n'+
            '  ! assign dimension to subroutine variables\n'+
            '  DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS, NTENS),'+
            '  DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),'+
            '  PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3), DROT(3, 3),'+
            '  DFGRD0(3, 3), DFGRD1(3, 3)\n\n')
    #endregion

    #region variable declaration
    f.write('! initialize algorithm variables\n'
            '  ! fixed variables: S_sigma, S_epsilon, beta, NB, NC, ND, NE, w_BA, w_CB, w_DC, w_ED, w_FE \n' +
            '  ! state dependent variables: sigma_NN_i1, B, C, D, E\n')
    f.write('  integer                                  :: i, j, k, l, NB, NC, n_dim\n'.format(n_dim) +
            '  real(DP)                                 :: S_sigma, S_epsilon, beta\n'
            '  real(DP), dimension ({})                 :: sigma\n'.format(n_dim) +
            '  real(DP), dimension ({}, 1)              :: sigma_NN_i1, sigma_out\n'.format(n_dim) +
            '  real(DP), dimension ({}, {})             :: Dmat, Dmatrx\n'.format(n_dim, n_dim) +
            '  real(DP), dimension ({}, 1)              :: B, b_BA\n'.format(w_BA.shape[0])+
            '  real(DP), dimension ({}, 1)              :: C, b_final\n'.format(b_final.shape[0])+
            '  real(DP), dimension ({}, 1)              :: network_input\n'.format(w_BA.shape[1]) +
            '  real(DP), dimension ({}, {})             :: w_BA\n'.format(w_BA.shape[0], w_BA.shape[1])+
            '  real(DP), dimension ({}, {})             :: w_final\n'.format(w_final.shape[0], w_final.shape[1]) +
            '  real(DP), dimension (:,:), allocatable   :: w_CB, b_CB, w_DC, b_DC, w_ED, b_ED\n'
            '  real(DP), dimension (:,:), allocatable   :: D, E, x, net_input, output\n\n')
    #endregion

    #region fixed variable assignment
    f.write('  ! fixed variable assignment\n')
    f.write('  CMNAME = "{}"\n  n_dim = 3\n'.format(matname))
    f.write('  S_sigma = {}\n'.format(normalisationParams.loc['max','Sigma_1']))
    f.write('  S_epsilon = {}\n'.format(normalisationParams.loc['max','Epsilon_1']))
    f.write('  S_dsigma = {}\n'.format(normalisationParams.loc['max','dSigma_1']))
    f.write('  S_depsilon = {}\n'.format(normalisationParams.loc['max','dEpsilon_1']))
    f.write('  beta = 1\n')
    f.write('  NB = {}\n'.format(model.layers[1].input_shape[1]))
    if W.shape[0] >= 3: f.write('  NC = {}\n'.format(model.layers[2].input_shape[1]))
    if W.shape[0] >= 4: f.write('  ND = {}\n'.format(model.layers[3].input_shape[1]))
    if W.shape[0] >= 5: f.write('  NE = {}\n'.format(model.layers[4].input_shape[1]))
    #endregion

    #region state dependent variable assignment
    f.write('  ! state dependent variable assignment\n'+
            '  do i = 1, n_dim\n'+
            '    network_input(i,1) = STRESS(i)/S_sigma\n'+
            '    network_input(i+n_dim,1) = STRAN(i)/S_epsilon\n'+
            '    network_input(i+n_dim*2,1) = DSTRAN(i)/S_depsilon\n'+
            '  end do\n\n')
    #endregion

    #region weight matrices
    for i in range(w_BA.shape[0]): f.write('  w_BA({},:) = {}\n'.format(i+1, ar2str(w_BA[i, :])))
    f.write('  b_BA(:,1) = {}\n'.format(ar2str(b_BA)))
    for i in range(w_final.shape[0]): f.write('  w_final({},:) = {}\n'.format(i+1, ar2str(w_final[i, :])))
    f.write('  b_final(:,1) = {}\n'.format(ar2str(b_final)))
    if W.shape[0] >= 3:
        w_CB, b_CB = W[1][0].T, W[1][1].reshape(-1,1)
        f.write('  allocate ( w_CB({},{}) )\n'.format(w_CB.shape[0], w_CB.shape[1]))
        for i in range(w_CB.shape[0]): f.write('  w_CB({},:) = {}\n'.format(i+1, ar2str(w_CB[i, :])))
        f.write('  allocate ( b_CB({},1) )\n'.format(w_CB.shape[0])+ '  b_CB(:, 1) = {}\n'.format(ar2str(b_CB)))
    if W.shape[0] >= 4:
        w_DC, b_DC = W[2][0].T, W[2][1].reshape(-1,1)
        f.write('  allocate ( w_DC({},{}) )\n'.format(w_DC.shape[0], w_DC.shape[1]))
        for i in range(w_DC.shape[0]): f.write('  w_DC({},:) = {}\n'.format(i+1, ar2str(w_DC[i, :])))
        f.write('  allocate ( b_DC({},1) )\n'.format(w_DC.shape[0])+ '  b_DC(:, 1) = {}\n'.format(ar2str(b_DC)))
    if W.shape[0] >= 5:
        w_ED, b_ED = W[3][0].T, W[3][1].reshape(-1, 1)
        f.write('  allocate ( w_ED({},{}) )\n'.format(w_ED.shape[0], w_ED.shape[1]))
        for i in range(w_ED.shape[0]): f.write('  w_ED({},:) = {}\n'.format(i+1, ar2str(w_ED[i, :])))
        f.write('  allocate ( b_ED({},1) )\n'.format(w_ED.shape[0])+ '  b_ED(:, 1) = {}\n'.format(ar2str(b_ED)))

    #endregion

    #region compute intermediate values of neural network
    f.write('\n\n  ! Neural network calculation\n'+
            '  ! First hidden layer calculation\n'
            '  allocate ( x ({}, 1))        ! neuron input\n'.format(w_BA.shape[1])+
            '  allocate ( net_input({}, 1)) ! shape same as w_BA.shape[0]\n'.format(w_BA.shape[0])+
            '  allocate ( output({}, 1))    ! shape same as w_BA.shape[0]\n'.format(w_BA.shape[0])+
            '  ! first layer calculation\n'+
            '  x = network_input\n'+
            '  net_input = matmul(w_BA,x) + b_BA\n'+
            '  output = max(net_input,ZERO)   ! ReLu function\n'+
            '  B = output\n')
    if W.shape[0] >= 3:
        f.write('  ! second hidden layer calculation\n'+
                '  deallocate ( net_input)\n'+
                '  deallocate (x)\n'+
                '  allocate ( x({}, 1))          ! shape same as w_CB.shape[1]\n'.format(w_CB.shape[1])+
                '  allocate ( net_input({}, 1))  ! shape same as w_CB.shape[0]\n'.format(w_CB.shape[0])+
                '  x = output\n'+
                '  deallocate (output)\n'+
                '  allocate ( output({}, 1))     ! shape same as w_CB.shape[0]\n'.format(w_CB.shape[0])+
                '  net_input = matmul(w_CB, x) + b_CB\n'+
                '  output = max(net_input,ZERO)    ! ReLu function\n'+
                '  C = output\n'+
                '  deallocate(w_CB)\n'+
                '  deallocate(b_CB)\n')
    if W.shape[0] >= 4:
        f.write('  ! third hidden layer calculation\n'+
                '  deallocate ( net_input)\n'+
                '  deallocate (x)\n'+
                '  allocate ( x({}, 1))          ! shape same as w_DC.shape[1]\n'.format(w_DC.shape[1])+
                '  allocate ( net_input({}, 1))  ! shape same as w_DC.shape[0]\n'.format(w_DC.shape[0])+
                '  x = output\n'+
                '  deallocate (output)\n'+
                '  allocate ( output({}, 1))     ! shape same as w_DC.shape[0]\n'.format(w_DC.shape[0])+
                '  net_input = matmul(w_DC, x) + b_DC\n'+
                '  output = max(net_input,ZERO)    ! ReLu function\n'+
                '  D = output\n'+
                '  deallocate(w_DC)\n'+
                '  deallocate(b_DC)\n')
    if W.shape[0] >= 5:
        f.write('  ! fourth hidden layer calculation\n' +
                '  deallocate ( net_input)\n' +
                '  deallocate (x)\n' +
                '  allocate ( x({}, 1))         ! shape same as w_ED.shape[1]\n'.format(w_ED.shape[1]) +
                '  allocate ( net_input({}, 1)) ! shape same as w_ED.shape[0]\n'.format(w_ED.shape[0]) +
                '  x = output\n' +
                '  deallocate (output)\n' +
                '  allocate ( output({}, 1))    ! shape same as w_ED.shape[0]\n'.format(w_ED.shape[0]) +
                '  net_input = matmul(w_ED, x) + b_ED\n' +
                '  output = max(net_input,ZERO)   ! ReLu function\n' +
                '  E = output\n'+
                '  deallocate(w_ED)\n'+
                '  deallocate(b_ED)\n')
    f.write('  ! output layer calculation\n' +
            '  deallocate (net_input)\n' +
            '  deallocate (x)\n' +
            '  allocate ( x({}, 1))             ! shape same as w_final.shape[0]\n'.format((W[-1][0].T).shape[1]) +
            '  allocate ( net_input({}, 1))     ! shape same as w_final.shape[1]\n'.format((W[-1][0].T).shape[0]) +
            '  x = output\n' +
            '  deallocate (output)\n' +
            '  allocate ( output({}, 1))        ! shape same as w_final.shape[1]\n'.format((W[-1][0].T).shape[0]) +
            '  net_input = matmul(w_final, x) + b_final\n' +
            '  output = tanh(net_input)\n' +
            '  sigma_NN_i1 = output\n'
            'deallocate(output)\n'+
            'deallocate(x)\n'+
            'deallocate(net_input)')
    #endregion

    #region compute stiffness matrix
    f.write('\n\n  ! Stiffness matrix calculation\n')

    if solver == 'hashash':
        f.write('  real(DP) :: Part0, PartA, PartA_A, PartA_B\n')
        # Specific routine for [relu, relu, tanh] network. Has to be re-derived for other structures
        if W.shape[0] == 3:
            f.write('  do i = 1, n_dim\n'+
                    '    do j = 1, n_dim\n'+
                    '      Part0 = S_sigma / S_epsilon * beta ** 3\n'+
                    '      PartAB = 0.\n'+
                    '      do k = 1, NC\n'+
                    '        PartA = ((1.-(sigma_NN_i1(i))**2)*w_DC(i, k))  ! represents the final tanh(...) layer\n'+
                    '        PartB = 0.\n'+
                    '        do l = 1, NB\n'+
                    '          PartB = PartB + min(ceiling(max(C[k],0.)),1.)*w_CB(k, l) * '
                                                'min(ceiling(max(B[k],0.)),1.)*w_BA(l, j) ! represents ReLu(...) layers\n'+
                    '          ! tanh version: PartA_B = ((1-(C(k))**2)*w_CB(l, k)) * ((1-(B(l))**2)*w_BA(l, j))\n'+
                    '        end do\n'+
                    '        PartAB = PartAB + PartA * PartB\n'+
                    '      end do\n'+
                    '      Dmat(i,j) = Part0 * PartAB\n'
                    '    end do\n'+
                    '  end do\n\n')
        elif W.shape[0] == 4:
            f.write('  do i = 1, n_dim\n' +
                    '    do j = 1, n_dim\n' +
                    '      Part0 = S_sigma / S_epsilon * beta\n' +
                    '      PartABC = 0.\n' +
                    '      do k = 1, ND\n' +
                    '        PartA = ((1.-(sigma_NN_i1(i))**2)*w_ED(k, i))  ! represents the final tanh(...) layer\n' +
                    '        PartB = 0.\n' +
                    '        PartBC = 0.\n'
                    '        do l = 1, NC\n' +
                    '          PartB = PartB + min(ceiling(max(D[k],0.)),1.)*w_DC(l, k) \n'+
                    '          PartC = 0.\n'+
                    '          do m = 1, NB\n'+
                    '            PartC = PartC + min(ceiling(max(C[l],0.)),1.)*w_CB(l, m) *'+
                                            'min(ceiling(max(B[m],0.)),1.)*w_BA(m, j) \n' +
                    '          end do\n'+
                    '          PartBC = PartBC + PartB * PartC\n'
                    '        end do\n' +
                    '        PartABC = PartABC + PartA * PartBC\n' +
                    '      end do\n' +
                    '      Dmat(i,j) = Part0 * PartABC\n'
                    '    end do\n' +
                    '  end do\n\n')
        else:
            raise ModuleNotFoundError('stiffness matrix computation not available for '+
                                      'network with {} layers'.format(W.shape[0]))
    elif solver == 'linear': # NOT DONE NOT DONE NOT DONE
        f.write('  real(DP), dimension(n_dim, n_dim) :: Amat, bmat\n'+
                '  real(DP), dimension(n_dim,n_dim) :: strain_state\n'+
                '  if n_dim>=2: \n'+
                '    strain_state(1,1) = network_input(n_dim+1,1)*4.\n'+
                '    strain_state(1,2) = network_input(n_dim+2,1)/2.\n'+
                '    strain_state(2,1) = network_input(n_dim+1,1)/2.\n'+
                '    strain_state(2,2) = network_input(n_dim+2,1)*4.\n'+
                '  else if n_dim>=3:\n'+
                '    strain_state(3,1) = network_input(n_dim+3,1)*4.\n'+
                '    strain_state(3,2) = network_input(n_dim+3,1)/2.\n'+
                '    strain_state(3,3) = network_input(n_dim+3,1)*2.\n'+
                '    strain_state(1,3) = network_input(n_dim+3,1)*2.\n'+
                '    strain_state(2,3) = network_input(n_dim+3,1)*2.\n'+
                '  end if\n'+
                '  \n'+
                '  \n'+
                '  \n'+
                '  \n')
    elif solver == 'secant':
        f.write('  ! secant stiffness matrix\n'+
                '  do i = 1, n_dim\n'+
                '    do j = 1, n_dim\n'+
                '      if (abs(stran(j)+dstran(j))<SMALLNUM) then\n'+
                '      Dmat(i, j) = stress(i) / SMALLNUM   ! 100._dp\n'+
                '      else\n'+
                '      Dmat(i, j) = stress(i) / stran(j)\n'+
                '      end if\n'+
                '    end do\n'+
                '  end do\n')
    else: raise ModuleNotFoundError('Solver type not implemented')
    #endregion

    #region backmatter
    f.write('  ! in the end, pass Umat and sigma to Abaqus UMAT Dmatrx and sigma_out\n'+
            '   DDSDDE   = Dmat\n'+
            '   STRESS   = STRESS + sigma_NN_i1(:,1)*S_dsigma\n\n')
    f.write('end subroutine umat\n')

    f.close()
    #endregion

    #region cap lines that are too long
    f = open(fname)
    lines = f.readlines()
    i = 0
    cap = 130
    while i != len(lines):
        if len(lines[i]) >= cap:
            part1, part2 = lines[i][:cap]+'&\n', '&'+lines[i][cap:]
            lines[i] = part1
            lines.insert(i+1, part2)
        i += 1
    f.close()
    f = open(fname, 'w')
    for line in lines:
        f.write(line)
    f.close()
    #endregion

    return None


#region miscellaneous stuff
# -------------------

# helper function to format plot
def plotStyle(ax, x_ticks, y_ticks, *args):
    if len(args) == 0 or 'standard' in args:
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

        plt.tick_params(axis="both", which="both", bottom="on", top="off", labelbottom="on", left="on", right="off",
                        labelleft="on")
        plt.axvline(c='k', linewidth=.5)
        plt.axhline(c='k', linewidth=.5)
        ax.set_xticks(x_ticks)
        ax.set_yticks(y_ticks)

        plt.xlabel('True strain [-]', fontsize=11)
        plt.ylabel('True stress [MPa]', fontsize=11)
    return None

# helper function for putting arrays into Fortran (free-form)
def ar2str(array):

    return '(/'+np.array2string(array.T, separator=',').replace(' ', '').replace('[', '').replace(']', '').\
        replace('\n', '')+'/)'

# helper function to visualise network structure
def plot_network(model):
    plot_model(model, to_file='model.png', show_shapes=True)
    return None

# helper function to call plot colors
def load_plotColors():
    tableau  = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
                (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
                (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
                (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
                (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

    for i in range(len(tableau)):
        r, g, b = tableau[i]
        tableau[i] = (r / 255., g / 255., b / 255.)
    return tableau

def ReLu(x): return x * (x > 0)

def tanh(x): return (np.exp(x)-np.exp(-x))/(np.exp(x)+np.exp(-x))

# -------------------
#endregion
# =================================
# Author: T.H.E. Gulikers
# Main file
# Description:
# Aluminium elementary brick, 2D, version 1. First iteration with purpose to demonstrate a parametrised model and
# data set generator with abaqus
# =================================

#region initialise modules
import sys
import os
import time
tick = time.time()

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from abaqusConstants import *
from abaqus import *
import regionToolset
import numpy as np
import odbAccess

# set viewport
session.viewports['Viewport: 1'].setValues(displayedObject = None)
#endregion

#region parameter inputs
path = 'C:/Users/tom/Documents/AQS_models/'
currentfilename = 'T_steelBrick_2_1_1_main.py'
materialChoice = 'SteelBoyang'

lengths  = [10]         # lengths of panel to be considered
widths   = [10]         # widths of panel to be considered
depths   = [1]          # depths of panel to be considered

minstress_2 = 450.      # [MPa] Minimum loading step in direction 1 to consider in creating the data set
maxstress_2 = 450.      # [MPa] Maximum loading step in direction 1 to consider in creating the data set
steps_2     = 1         # amount of loading points in direction 1 to consider in creating the data set

minstress_1 = 1e-10     # [MPa] Minimum loading step in direction 2 to consider in creating the data set
maxstress_1 = 450.      # [MPa] Maximum loading step in direction 2 to consider in creating the data set
steps_1     = 7         # amount of loading points in direction 2 to consider in creating the data set

dsigmas    = [1., 2., 4., 8., 12.]    # [MPa] increment to be used in Abaqus loading step. set 0 for default of 1 MPa
seedSizeGlobal = 2.
#endregion

#region process parameters and delete leftover models
modelGenerator, resultsExtractor = currentfilename.replace('main', 'modelGenerator'), currentfilename.replace('main', 'resultsExtractor')
sigma_1_datalist = np.linspace(minstress_1, maxstress_1, steps_1)       # generate list of stresses in direction 1
sigma_2_datalist = np.linspace(minstress_2, maxstress_2, steps_2)       # generate list of stresses in direction 2

mdb.Model(name='Model-1')
for modelName in mdb.models.keys():
    if modelName != 'Model-1':
        del mdb.models[modelName]
for jobname in mdb.jobs.keys(): del mdb.jobs[jobname]

print "ready to initialise models. T: ", time.time()-tick
#endregion

#region create models
for length in lengths:
    for width in widths:
        for depth in depths:
            for sigma_1 in sigma_1_datalist:
                for sigma_2 in sigma_2_datalist:
                    for dsigma in dsigmas:
                        execfile(path + modelGenerator)

toc = time.time()
print len(mdb.models.keys()), " Models created. T: ", toc-tick
#endregion

if True:
    #region run jobs
    for jobname in mdb.jobs.keys():
        mdb.jobs[jobname].submit()
        mdb.jobs[jobname].waitForCompletion()
        print "Finished job ", jobname, " T: ", time.time() - toc
        toc = time.time()
    #endregion

    #region extract data from odb's
    for jobname in mdb.jobs.keys():
        execfile(path+resultsExtractor)
    #endregion

print "Abaqus has finished. Total runtime: ", time.time()-tick
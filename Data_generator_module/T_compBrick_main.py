# =================================
# Author: T.H.E. Gulikers
# Main file
# Description:
# Composite plate with hole, 2D, version 2, iteration 3.
# purpose: Different way of generating input parameters in order to compare to T_compBrick_2_2_2
# =================================

#region initialise modules
import sys
import os
import time
import numpy as np
from itertools import product

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
import odbAccess

os.environ["ABA_ACCELERATOR_TYPE"]="PLATFORM_CUDA" # Nvidia GPU configuration

tick = time.time()                                                      # reference time

session.viewports['Viewport: 1'].setValues(displayedObject = None)      # set viewport
#endregion

#region parameter inputs
path = 'C:/Users/tom/Documents/AQS_models/'
currentfilename = 'T_compBrick_2_2_3_main.py'
materialChoice = 'Composite_boyang1'

lengths  = (10,)                 # [mm] lengths of panel to be considered, y-direction
widths   = (10,)                 # [mm] widths of panel to be considered, x-direction
depths   = (1,)                  # [mm] depths of panel to be considered
position_ellipse = (5., 5.)     # [mm] width, length indicating ellipse centre
dimension_ellipse = (1., 2.)    # [mm] half of the hole size axis along [x, y] axis.

layup = (45., 90., -45., 0.)    # [deg] ply specification
makeSymmetric = True            # creates symmetric layup if True


# provide a minimum and maximum strain and amount of steps
minstrain_1 = 0.0e-2        # [-] Minimum loading step in direction 1 to consider in creating the data set
maxstrain_1 = 2.0e-2        # [-] Maximum loading step in direction 1 to consider in creating the data set
steps_1     = 4             # amount of loading points in direction 1 to consider in creating the data set

minstrain_2 = 0.0e-2        # [-] Minimum loading step in direction 2 to consider in creating the data set
maxstrain_2 = 2.0e-2        # [-] Maximum loading step in direction 2 to consider in creating the data set
steps_2     = 4             # amount of loading points in direction 2 to consider in creating the data set

minstrain_12 = 0.0e-2       # [-] Minimum loading step in direction 2 to consider in creating the data set
maxstrain_12 = 2.0e-2       # [-] Maximum loading step in direction 2 to consider in creating the data set
steps_12     = 4            # amount of loading points in direction 2 to consider in creating the data set

# provide a minimum increment size and list of multipliers (simulation step time is 1)
dt = 5e-3                   # [-] increment to be used in Abaqus loading step. 0.0001
dt_multipliers = (1,3,8) # [-] tuple of integers. used to extract results with different step sizes

seedsOnEdge = 15            # amount of seeds per outer edge. The elliptic hole edge has 4x this number of seeds
#endregion


#region process parameters and delete previous models
modelGenerator, resultsExtractor = currentfilename.replace('main', 'modelGenerator'), currentfilename.replace('main', 'resultsExtractor')
epsilon_1_datalist = np.linspace(minstrain_1, maxstrain_1, steps_1)       # generate list of strains in direction 1
epsilon_2_datalist = np.linspace(minstrain_2, maxstrain_2, steps_2)       # generate list of strains in direction 2
epsilon_12_datalist = np.linspace(minstrain_12, maxstrain_12, steps_12)       # generate list of strains in direction 12

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
            for epsilon_1 in epsilon_1_datalist:
                for epsilon_2 in epsilon_2_datalist:
                    for epsilon_12 in epsilon_12_datalist:
                        if epsilon_1 != maxstrain_1 and epsilon_2 != maxstrain_2 and epsilon_12 != maxstrain_12:
                            continue
                        else: execfile(path + modelGenerator)

toc = time.time()
print len(mdb.models.keys()), " Models created. T: ", toc-tick
#endregion


# create jobs
if False:
    for jobname in mdb.jobs.keys():
        mdb.jobs[jobname].submit()
        mdb.jobs[jobname].waitForCompletion()
        print "Finished ", jobname, " T: ", time.time() - toc
        toc = time.time()

# extract data from odb's
if False:
    for jobname in mdb.jobs.keys():
       execfile(path+resultsExtractor)


print "Abaqus has finished. Total runtime: ", time.time()-tick

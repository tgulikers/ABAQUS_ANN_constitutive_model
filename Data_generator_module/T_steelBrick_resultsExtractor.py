# =================================
# Author: T.H.E. Gulikers
# Results analysis file
# Description:
# Aluminium elementary brick, 2D, version 1. First iteration with purpose to demonstrate a parametrised model and
# data set generator with abaqus
# =================================

# open .odb object
odbName = os.getcwd()+'\\'+jobname+'.odb'
odbObject = session.openOdb(name=odbName)


#region analyse dimensions through the jobname
locL, locW, locD = jobname.find('L'), jobname.find('W'), jobname.find('D')      # find key letters indicating length, width, depth
dimLst, j = ['', '', ''], 0
for loc in [locL, locW, locD]:
    current, i, num = jobname[loc+1], 2, []     # initialise loop parameters. current (currently considered string entry), i (index of current), num (list of valid numbers)
    while current != '_':           # find how many numbers there are after each letter
        num.append(current)
        current = jobname[loc+i]
        i += 1

    for k in num:
        dimLst[j] = int(str(dimLst[j])+k)
    j += 1
length, width, depth = [float(i) for i in dimLst]
#endregion

#region plot deformed state
viewport_deformed = session.viewports['Viewport: 1']
viewport_deformed.setValues(displayedObject=odbObject)
viewport_deformed.odbDisplay.display.setValues(plotState=(UNDEFORMED, DEFORMED, ))
viewport_deformed.odbDisplay.commonOptions.setValues(nodeLabels=ON)
viewport_deformed.odbDisplay.commonOptions.setValues(elemLabels=ON)
viewport_deformed.setValues(origin=(0.0, 0.0))
#endregion

#region label relevant node sets
setBottom   = odbObject.rootAssembly.nodeSets['EDGE_BOTTOM']
setLeft     = odbObject.rootAssembly.nodeSets['EDGE_LEFT']
setTop      = odbObject.rootAssembly.nodeSets['EDGE_TOP']
setRight    = odbObject.rootAssembly.nodeSets['EDGE_RIGHT']
#endregion

#region Extract field outputs and write to output file
reportFileName = os.getcwd()+'\\'+jobname.strip('Job_')+'.txt'
try:
    os.remove(reportFileName)
except:
    pass

reportFile = open(reportFileName, 'w')
reportTitle = 'Sigma_1,Sigma_2,Epsilon_1,Epsilon_2,dSigma_1,dSigma_2,dEpsilon_1,dEpsilon_2'

step = odbObject.steps.values()[0]
dataArr = np.array([[0, 0, 0, 0, 0, 0, 0, 0]])
sigma_1prev, sigma_2prev = 0., 0.
eps_1prev, eps_2prev = 0., 0.
for frame in step.frames:
    currentTimeValue = frame.frameValue     # current time value
    RFfield = frame.fieldOutputs['RF']
    Ufield = frame.fieldOutputs['U']

    RF_b, RF_l, U_t, U_r = RFfield.getSubset(region=setBottom), RFfield.getSubset(region=setLeft),\
                                Ufield.getSubset(region=setTop), Ufield.getSubset(region=setRight)

    RFval_b, RFval_l, Uval_t, Uval_r = RF_b.values, RF_l.values, U_t.values, U_r.values
    numNodes_b, numNodes_l, numNodes_t, numNodes_r = len(RFval_b), len(RFval_l), len(Uval_t), len(Uval_r)

    # sum edge forces and displacements
    RF2_b, RF1_l, U2_t, U1_r = 0., 0., 0., 0.
    for nodeNum in range(numNodes_b): RF2_b += RFval_b[nodeNum].data[1]
    for nodeNum in range(numNodes_l): RF1_l += RFval_l[nodeNum].data[0]
    for nodeNum in range(numNodes_t): U2_t  += Uval_t[nodeNum].data[1]
    for nodeNum in range(numNodes_r): U1_r  += Uval_r[nodeNum].data[0]

    # compute current stresses and strains
    eps_1cur, eps_2cur      = np.log(1.+U1_r/width/numNodes_r), np.log(1.+U2_t/length/numNodes_t)         # true strain
    sigma_1cur, sigma_2cur  = -RF1_l/width/depth*(1.+eps_1cur), -RF2_b/length/depth *(1.+eps_2cur)        # true stress

    # compute stress and strain increments
    dsigma_1, dsigma_2      = sigma_1cur - sigma_1prev, sigma_2cur-sigma_2prev
    deps_1, deps_2          = eps_1cur-eps_1prev, eps_2cur-eps_2prev

    # database entry
    arrEntry = np.array([[sigma_1prev, sigma_2prev, eps_1prev, eps_2prev, dsigma_1, dsigma_2, deps_1, deps_2]])
    dataArr = np.concatenate( (dataArr, arrEntry ), axis=0)

    # update values
    sigma_1prev, sigma_2prev, eps_1prev, eps_2prev = sigma_1cur, sigma_2cur, eps_1cur, eps_2cur

dataArr = dataArr[2:]       # remove all zero rows
reportFile.write(reportTitle+'\n')
np.savetxt(reportFile, dataArr, delimiter = ',')
reportFile.close()
#endregion

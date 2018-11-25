# =================================
# Author: T.H.E. Gulikers
# Results analysis file
# Description:
# Composite plate with hole, 2D, version 2, iteration 3.
# purpose: Changed how strain is extracted
# =================================

# open .odb object
odbName = os.getcwd()+'\\'+jobname+'.odb'
odbObject = session.openOdb(name=odbName)
odbInstance = odbObject.rootAssembly.instances['Part Instance']

#region analyse dimensions through the jobname
locL, locW, locD = jobname.find('L'), jobname.find('W'), jobname.find('D')      # find key letters indicating length, width, depth
dimLst, j = ['', '', ''], 0
for loc in [locL, locW, locD]:
    current, i, num = jobname[loc+1], 2, []   # initialise loop parameters. current (currently considered string entry),
                                              #  i (index of current string entry), num (list of valid numbers)
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
# viewport_deformed.odbDisplay.commonOptions.setValues(nodeLabels=ON)
# viewport_deformed.odbDisplay.commonOptions.setValues(elemLabels=ON)
# viewport_deformed.setValues(origin=(0.0, 0.0))
#endregion

#region label relevant node sets
setBottom   = odbObject.rootAssembly.nodeSets['EDGE_BOTTOM']
setLeft     = odbObject.rootAssembly.nodeSets['EDGE_LEFT']
setTop      = odbObject.rootAssembly.nodeSets['EDGE_TOP']
setRight    = odbObject.rootAssembly.nodeSets['EDGE_RIGHT']
setNormal   = odbObject.rootAssembly.nodeSets['RP_NORMAL']
setShear    = odbObject.rootAssembly.nodeSets['RP_SHEAR']

# find corner nodes
nodeBottom_L = [i for i in range(len(setBottom.nodes[0])) if setBottom.nodes[0][i].coordinates[0] == 0.][0]
nodeBottom_R = [i for i in range(len(setBottom.nodes[0])) if setBottom.nodes[0][i].coordinates[0] == width][0]
nodeLeft_B   = [i for i in range(len(setLeft.nodes[0])) if setLeft.nodes[0][i].coordinates[1] == 0.][0]
nodeLeft_T   = [i for i in range(len(setLeft.nodes[0])) if setLeft.nodes[0][i].coordinates[1] == length][0]
#endregion


# create multiple instances of simulation with different increments
for multiplier in dt_multipliers:
    dt_current = dt*multiplier                                  # increment size
    k = 0                                                       # counter

    # region Extract field outputs
    step = odbObject.steps.values()[0]                          # list of step 1 frames
    sig_1prev, sig_2prev, sig_12prev = 0., 0., 0.               # initiate variables to store frame stress
    eps_1prev, eps_2prev, eps_12prev = 0., 0., 0.               # initiate variables to store frame strain

    dataArr = np.zeros((int(1./dt_current)+1, 12))
    for frame in step.frames:
        currentTimeValue = frame.frameValue     # current time value

        if currentTimeValue >= dt_current*k:    # this ensures a minimum increment size in the output file
            pass
        else:
            continue

        # store reaction-force and displacement field
        RFfield = frame.fieldOutputs['RF']
        Ufield  = frame.fieldOutputs['U']

        # extract lists of forces and displacements at edges
        RF_1, RF_2, RF_12       = RFfield.getSubset(region=setNormal).values[0].data[0],\
                                 RFfield.getSubset(region=setNormal).values[0].data[1],\
                                 [RFfield.getSubset(region=setShear).values[0].data[0],
                                 RFfield.getSubset(region=setShear).values[0].data[1]]
        U_1, U_2, U_12          = Ufield.getSubset(region=setNormal).values[0].data[0], \
                                    Ufield.getSubset(region=setNormal).values[0].data[1], \
                                    [Ufield.getSubset(region=setShear).values[0].data[0],
                                    Ufield.getSubset(region=setShear).values[0].data[1]]


        # compute current stresses and strains and convert engineering to true values
        eps_1cur  = np.log(1.+(U_1/width))
        eps_2cur  = np.log(1.+(U_2/length))
        eps_12cur = U_12[0]/length

        sig_1cur  = RF_1/width/depth * (np.exp(eps_1cur))
        sig_2cur  = RF_2/length/depth * (np.exp(eps_2cur))
        sig_12cur = RF_12[0]/width/depth+RF_12[1]/length/depth


        # compute stress and strain increments
        dsig_1, dsig_2, dsig_12 = sig_1cur - sig_1prev, sig_2cur - sig_2prev, sig_12cur - sig_12prev
        deps_1, deps_2, deps_12 = eps_1cur - eps_1prev, eps_2cur - eps_2prev, eps_12cur - eps_12prev

        # database entry
        dataArr[k, :] = np.array([[sig_1prev, sig_2prev, sig_12prev, eps_1prev, eps_2prev, eps_12prev, dsig_1, dsig_2,
                                   dsig_12, deps_1, deps_2, deps_12]])

        # update values
        sig_1prev, sig_2prev, sig_12prev = sig_1cur, sig_2cur, sig_12cur
        eps_1prev, eps_2prev, eps_12prev = eps_1cur, eps_2cur, eps_12cur
        k += 1

    dataArr = dataArr[1:]       # remove all zero rows

    # write to file
    reportFileName = os.getcwd()+'\\'+jobname.strip('Job_').strip('1')+str(int(dt*1e3*multiplier))+'.txt'
    reportFileName = reportFileName.replace('T_compBrick_2_2_3','T_compBrick')
    reportFile = open(reportFileName, 'w')
    reportTitle = 'Sigma_1,Sigma_2,Sigma_12,Epsilon_1,Epsilon_2,Epsilon_12,dSigma_1,dSigma_2,dSigma_12,dEpsilon_1,' \
                  'dEpsilon_2,dEpsilon_12'
    reportFile.write(reportTitle+'\n')
    np.savetxt(reportFile, dataArr, delimiter = ',')
    reportFile.close()
#endregion
# =================================
# Author: T.H.E. Gulikers
# Model generator
# Input: dimensions and loading
# Process: defines part, sections, materials, mesh, assembly, steps, constraints, loads and job
# Output: Model and job file (not submitted)
# Description:
# Composite plate with hole, 2D, version 2, iteration 3.
# purpose: Fix application of strain
# =================================

#region create the model
modelName = currentfilename.strip('_2_2_2_main.py') + '_L' + str(length) + '_W' +  str(width) + '_D' +  str(depth) + \
            '_Ea' +  str(int(epsilon_1*1e6)) + '_Eb' + str(int(epsilon_2*1e6)) + '_Es' + str(int(epsilon_12*1e6)) + '_dE' + str(1)

if 'Model-1' in mdb.models.keys():
    mdb.models.changeKey(fromName='Model-1', toName=modelName)
else:
    mdb.Model(name=modelName)

currentModel = mdb.models[modelName]
#endregion


#region create part
currentSketch = currentModel.ConstrainedSketch(name='2D_sketch', sheetSize=float(np.max((width, length, depth))))
currentSketch.rectangle((0., 0.), (width, length))

# draw hole
currentSketch.EllipseByCenterPerimeter(center= (width/2., length/2.), axisPoint1 =(width/2., length/2.-dimension_ellipse[1]),
                                                    axisPoint2 = (width/2.-dimension_ellipse[0], length/2.))

# part from sketch
currentPart = currentModel.Part(name='part_3D', dimensionality=THREE_D, type=DEFORMABLE_BODY)
currentPart.BaseShell(sketch=currentSketch)
del currentModel.sketches['2D_sketch']

# coordinate system
CSYS_base = currentPart.DatumCsysByThreePoints(coordSysType= CARTESIAN, line1=(1.0, 0.0, 0.0),
                                            line2=(0.0, 1.0, 0.0), name='csys_base', origin=(0.0, 0.0, 0.0))
#endregion


#region import materials from database
allMats = False
execfile(path+'MatDatabase.py')
#endregion


#region create composite section
sectionRegion = (currentPart.faces.findAt((width*0.99, length*0.99, 0.),),) # location on the plate

# composite layup object initiation
currentSection = currentPart.CompositeLayup(description='', elementType=SHELL, name='QI_section',
            offsetType=MIDDLE_SURFACE, symmetric=makeSymmetric, thicknessAssignment=FROM_SECTION)
currentSection.Section(integrationRule=SIMPSON, poissonDefinition=DEFAULT, preIntegrate=OFF, thicknessType=UNIFORM,
                       useDensity=OFF)
currentSection.ReferenceOrientation(additionalRotationType=ROTATION_NONE, angle=0.0, axis=AXIS_3, fieldName='',
    localCsys=None, orientationType=GLOBAL)

# create plies

# get ply thickness from plate depth and layup
if makeSymmetric: t_ply = depth/2./len(layup)
else: t_ply = depth/len(layup)
i = 0                                           # enumerate ply names
for angle_ply in layup:
    i += 1
    currentSection.CompositePly(thickness=t_ply, angle=angle_ply, axis=AXIS_3, material= materialChoice,
                                orientationType=CSYS, orientation=currentPart.datums[CSYS_base.id],
                                plyName='Ply-'+str(i), region=sectionRegion, thicknessType=SPECIFY_THICKNESS)
#endregion


#region create assembly and define edges
currentAssembly = currentModel.rootAssembly
currentAssembly.Instance(name='Part Instance', part=currentPart, dependent=ON)
currentInstance = currentAssembly.instances['Part Instance']

# bottom, left, top, right, edge objects
edgefind = currentInstance.edges.findAt
edges_b = edgefind((( width/2., 0., 0.),),)
edges_l = edgefind(((0., length/2., 0.),),)
edges_t = edgefind(((width/2., length, 0.),),)
edges_r = edgefind(((width, length/2., 0.),),)
edges_ellipse = edgefind(((width/2.-dimension_ellipse[0], length/2.,0.),),)

# define sets from edges
currentAssembly.Set(name='edge_bottom', edges=edges_b)
currentAssembly.Set(name='edge_left',   edges=edges_l)
currentAssembly.Set(name='edge_top',    edges=edges_t)
currentAssembly.Set(name='edge_right',  edges=edges_r)
currentAssembly.Set(name='edge_circle',  edges=edges_ellipse)
#endregion


#region generate mesh (plane stress element)
currentMeshRegion = sectionRegion
meshElemType = ElemType(elemCode=S4, elemLibrary=STANDARD)

currentPart.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=seedsOnEdge)
currentPart.seedEdgeByNumber(constraint=FINER, edges=edges_ellipse, number=seedsOnEdge*4)
for edge in (edges_l, edges_r):  currentPart.seedEdgeByNumber(constraint=FINER, edges=edge, number=seedsOnEdge)
for edge in (edges_b, edges_t):  currentPart.seedEdgeByNumber(constraint=FIXED, edges=edge, number=seedsOnEdge)

currentPart.setMeshControls(algorithm=MEDIAL_AXIS, elemShape=QUAD, regions=currentMeshRegion)
currentPart.setElementType(regions=currentMeshRegion, elemTypes=(meshElemType,))
currentPart.generateMesh()
currentAssembly.regenerate()
#endregion


#region create step and field/history outputs
incFrac = dt
currentModel.StaticStep(name='ApplyLoad1', previous='Initial', description='Load is applied in this step',
                        initialInc=incFrac, maxInc=incFrac, maxNumInc=100000, nlgeom = ON, minInc=1e-15)
currentModel.steps['ApplyLoad1'].setValues(extrapolation=NONE, adaptiveDampingRatio=5e-04, continueDampingFactors=False,
                                           stabilizationMagnitude=5e-05, stabilizationMethod=DISSIPATED_ENERGY_FRACTION)
currentModel.steps['ApplyLoad1'].control.setValues(allowPropagation=OFF, resetDefaultValues=OFF,
                                                   timeIncrementation=(4.0, 8.0, 9.0, 16.0, 10.0, 4.0,
                                                                       12.0, 15.0, 6.0, 3.0, 50.0))

currentModel.fieldOutputRequests.changeKey(fromName='F-Output-1', toName='FieldOutputGlobal')
currentModel.fieldOutputRequests['FieldOutputGlobal'].setValues(variables=('S', 'E', 'PE', 'U', 'RF', 'CF'))
currentModel.FieldOutputRequest(createStepName='ApplyLoad1', layupLocationMethod=SPECIFIED, layupNames=(
    'Part Instance.QI_section', ), name='FieldOutputComposite', outputAtPlyBottom=False,
    outputAtPlyMid=True, outputAtPlyTop=False, rebar=EXCLUDE, variables=('S','MISES', 'E', 'PE', 'PEEQ', 'NE',
    'DAMAGEC', 'DAMAGET', 'DAMAGEFT','DAMAGEFC', 'DAMAGEMT', 'DAMAGEMC', 'DAMAGESHR', 'HSNFTCRT', 'HSNFCCRT',
    'HSNMTCRT', 'HSNMCCRT'))

currentModel.historyOutputRequests.changeKey(fromName='H-Output-1', toName='HistoryOutputs')
currentModel.historyOutputRequests['HistoryOutputs'].setValues(variables=PRESELECT)
#endregion


#region apply (periodic) boundary conditions

# find matching node pairs on model boundaries
pairs1, pairs2 = [], []                           # matching nodes for PBC in DoF 1 and DoF 2 (global coordinate system)
nodelst_dtype = [('label', 'i4'), ('x', 'f4'), ('y', 'f4'), ('z', 'f4')]  # data type specification needed for sorting
edgelabels = []

for set_name in ['edge_left', 'edge_right']:     # loop over direction 1
    for i_node in range(len(currentAssembly.sets[set_name].nodes)):         # loop over nodes in set
        label = currentAssembly.sets[set_name].nodes[i_node].label      # save node number and coordinates in list
        coords = currentAssembly.sets[set_name].nodes[i_node].coordinates
        pairs1.append([label]+list(coords))

        # define set (required for equation constraint)
        currentAssembly.SetFromNodeLabels(name='Node_' + str(label), nodeLabels=(('Part Instance', (label,)),))
        edgelabels.append(label)

for set_name in ['edge_bottom', 'edge_top']:     # loop over direction 2
    for i_node in range(len(currentAssembly.sets[set_name].nodes)):         # loop over nodes in set
        label = currentAssembly.sets[set_name].nodes[i_node].label      # save node number and coordinates in list
        coords = currentAssembly.sets[set_name].nodes[i_node].coordinates
        pairs2.append([label]+list(coords))

        # define set (required for equation constraint)
        currentAssembly.SetFromNodeLabels(name='Node_' + str(label), nodeLabels=(('Part Instance', (label,)),))
        edgelabels.append(label)

# convert lists of pairs to numpy arrays
pairs1, pairs2 = np.array([tuple(i) for i in pairs1], dtype=nodelst_dtype), np.array([tuple(i) for i in pairs2], dtype=nodelst_dtype)

# sort pair-lists based on coordinates such that the first half of the list is edge A and the second half is edge B
pairs1.sort(axis=0, order = ['x', 'y'])     # first x to divide the sides, then y to sort sequence on edge
pairs2.sort(axis=0, order = ['y', 'x'])     # first y to divide the sides, then x to sort sequence on edge
mp1, mp2 = len(pairs1)/2, len(pairs2)/2     # indices of mid-points, which is the division between edges

# define reference points to be used for loading and define a set for each
RP_normal_node = currentPart.Node(coordinates=(width*1.1, length/2., 0.)).setValues(label=90001)
RP_shear_node = currentPart.Node(coordinates=(width*-0.1, length/2., 0.)).setValues(label=90002)
currentAssembly.SetFromNodeLabels(name='RP_normal', nodeLabels=(('Part Instance', (90001,)),))
currentAssembly.SetFromNodeLabels(name='RP_shear', nodeLabels=(('Part Instance', (90002,)),))

# equation constraints
for i_pair in range(0, int(np.shape(pairs1)[0])/2):      # constraints on pairs in x-direction
    equationName = 'eq_'+str(pairs1[i_pair][0])+'_'+str(pairs1[i_pair+mp1][0])   # combination of node numbers

    # name, sequence of (Float, String, Int): (coefficient, Set name, DoF)
    currentModel.Equation('S'+equationName, ((1, 'Node_' + str(pairs1[i_pair][0]), 2),
			(-1, 'Node_' + str(pairs1[i_pair+mp1][0]), 2), (1, 'RP_shear', 1)))   # tangential

    if i_pair != 0:# and i_pair != int(np.shape(pairs1)[0]-2):
        currentModel.Equation('N'+equationName, ((1, 'Node_' + str(pairs1[i_pair][0]), 1),
			(-1, 'Node_' + str(pairs1[i_pair+mp1][0]), 1), (1, 'RP_normal', 1)))  # normal

for i_pair in range(0, int(np.shape(pairs2)[0])/2):      # constraints on pairs in y-direction
    equationName = 'eq_'+str(pairs2[i_pair][0])+'_'+str(pairs2[i_pair+mp2][0])   # combination of node numbers

    currentModel.Equation('S'+equationName, ((1, 'Node_' + str(pairs2[i_pair][0]), 1),
                          (-1, 'Node_' + str(pairs2[i_pair+mp2][0]), 1), (1, 'RP_shear', 2)))  # tangential

    # name, sequence of (Float, String, Int): (coefficient, Set name, DoF)
    if i_pair != 0:# and i_pair != int(np.shape(pairs2)[0]-2):
        currentModel.Equation('N'+equationName, ((1, 'Node_' + str(pairs2[i_pair][0]), 2),
			    (-1, 'Node_' + str(pairs2[i_pair+mp2][0]), 2), (1, 'RP_normal', 2)))  # normal
#endregion


#region apply load

# create meshNodeObjects from reference nodes (required for Region function)
meshNodeObjNormal = currentInstance.nodes.sequenceFromLabels((90001,))
meshNodeObjShear = currentInstance.nodes.sequenceFromLabels((90002,))

# apply normal loads
currentModel.DisplacementBC('other', createStepName='Initial',
                            region=Region(faces=currentInstance.faces.getSequenceFromMask(mask=('[#1 ]', ), )), u3=0.)

currentModel.DisplacementBC('load_RP_normal', createStepName='Initial', region=Region(nodes=meshNodeObjNormal), u1=0.,
                            u2=0.)

currentModel.boundaryConditions['load_RP_normal'].setValuesInStep('ApplyLoad1', u1 = epsilon_1*width,
                                                                  u2 =epsilon_2*length)

# apply shear load
currentModel.DisplacementBC('load_RP_shear', createStepName='Initial', region=Region(nodes=meshNodeObjShear), u1=0.,
                            u2=0.)

currentModel.boundaryConditions['load_RP_shear'].setValuesInStep('ApplyLoad1', u1 = epsilon_12*width,
                                                                 u2 = epsilon_12*length)

# find a node that is not on the edge and pin it such that the model doesn't fly away
for i in range(1,100000):
    if i not in edgelabels:
        break
currentModel.DisplacementBC('pin', createStepName='Initial', u1 = 0., u2 = 0.,
                             region=Region(nodes=currentInstance.nodes.sequenceFromLabels((i,))))
#endregion


#region create job
mdb.Job(name='Job_'+modelName, model=modelName, type=ANALYSIS, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE,
        description = 'Analysis of '+modelName, parallelizationMethodExplicit=DOMAIN, multiprocessingMode=DEFAULT,
        userSubroutine='', numCpus=5, numDomains=5, numGPUs=0, memory=70, memoryUnits=PERCENTAGE, echoPrint=OFF, modelPrint=OFF,
        contactPrint=OFF, historyPrint=OFF)
#endregion

# execfile('../AQS_models/T_compBrick_2_2_2_main.py')

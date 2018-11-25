# =================================
# Author: T.H.E. Gulikers
# Model generator
# Input: dimensions and loading
# Process: defines part, sections, materials, mesh, assembly, steps, constraints, loads and job
# Output: Model and job file (not submitted)
# Description:
# Aluminium elementary brick, 2D, version 1. First iteration with purpose to demonstrate a parametrised model and
# data set generator with abaqus
# =================================

#region create the model
modelName = currentfilename.strip('_2_1_1_main.py') + '_L' + str(length) + '_W' +  str(width) + '_D' +  str(depth) + \
            '_Sa' +  str(int(sigma_1)) + '_Sb' + str(int(sigma_2)) + '_dS' + str(int(dsigma))

if 'Model-1' in mdb.models.keys():
    mdb.models.changeKey(fromName='Model-1', toName=modelName)
else:
    mdb.Model(name=modelName)

currentModel = mdb.models[modelName]
#endregion


#region create part
currentSketch = currentModel.ConstrainedSketch(name='2D_sketch', sheetSize=float(np.max((length,width,depth))))
currentSketch.rectangle((length, 0),(0, width))

currentPart = currentModel.Part(name='part_2D', dimensionality=TWO_D_PLANAR, type=DEFORMABLE_BODY)
currentPart.BaseShell(sketch=currentSketch)
del currentModel.sketches['2D_sketch']
#endregion


#region import materials from database
allMats = False
execfile(path+'MatDatabase.py')
#endregion


#region create section
currentSection = currentModel.HomogeneousSolidSection(material=materialChoice, name='section_2D', thickness=depth)
sectionRegion = (currentPart.faces.findAt((length*0.99, width*0.99, 0.),),)
currentPart.SectionAssignment(region=sectionRegion, sectionName='section_2D', thicknessAssignment=FROM_SECTION)
#endregion


#region create assembly
currentAssembly = currentModel.rootAssembly
currentAssembly.Instance(name='Part Instance', part=currentPart, dependent=ON)
currentInstance = currentAssembly.instances['Part Instance']
#endregion


#region create step and field/history outputs
if dsigma==0: dsigma = 1.0
incFrac = round(1./np.max((sigma_1, sigma_2))*float(dsigma), 4)
currentModel.StaticStep(name='ApplyLoad1', previous='Initial', description='Load is applied in this step',
                        initialInc=incFrac, maxInc=incFrac, maxNumInc=5000, nlgeom = ON)

currentModel.fieldOutputRequests.changeKey(fromName='F-Output-1', toName='FieldOutputs')
currentModel.fieldOutputRequests['FieldOutputs'].setValues(variables=('S', 'E', 'PE', 'U', 'RF', 'CF'))

currentModel.historyOutputRequests.changeKey(fromName='H-Output-1', toName='HistoryOutputs')
currentModel.historyOutputRequests['HistoryOutputs'].setValues(variables=PRESELECT)
#endregion


#region apply load

# bottom, left, top, right edge objects
edgefind = currentInstance.edges.findAt
edges_b = edgefind((( width/2., 0., 0.),),)
edges_l = edgefind(((0., length/2., 0.),),)
edges_t = edgefind(((width/2., length, 0.),),)
edges_r = edgefind(((width, length/2., 0.),),)

# define sets from edges
currentAssembly.Set(name='edge_bottom', edges=edges_b)
currentAssembly.Set(name='edge_left',   edges=edges_l)
currentAssembly.Set(name='edge_top',    edges=edges_t)
currentAssembly.Set(name='edge_right',  edges=edges_r)

currentModel.Pressure(amplitude=UNSET, createStepName='ApplyLoad1', name='LoadRight', region=Region(side1Edges=edges_r,),
                      magnitude=-1.*sigma_1, distributionType=UNIFORM)
currentModel.Pressure(amplitude=UNSET, createStepName='ApplyLoad1', name='LoadTop', region=Region(side1Edges=edges_t,),
                      magnitude=-1.*sigma_2, distributionType=UNIFORM)

#endregion


#region apply boundary conditions
# pinned boundary conditions bottom and left
currentModel.YsymmBC(createStepName='Initial', name='YSymmBottom', region=Region(edges=edges_b,))
currentModel.XsymmBC(createStepName='Initial', name='XSymmLeft', region=Region(edges=edges_l,))
#endregion


#region generate mesh (plane stress element)
currentMeshRegion = sectionRegion
currentPart.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=seedSizeGlobal)
meshElemType = ElemType(elemCode=CPE4, elemLibrary=STANDARD)
currentPart.setElementType(regions=currentMeshRegion, elemTypes=(meshElemType,))
currentPart.generateMesh()
#endregion

#region create job
mdb.Job(name='Job_'+modelName, model=modelName, type=ANALYSIS, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE,
        description = 'Analysis of '+modelName, parallelizationMethodExplicit=DOMAIN, multiprocessingMode=DEFAULT,
        userSubroutine='', numCpus=7, numDomains=7, numGPUs=1, memory=70, memoryUnits=PERCENTAGE, echoPrint=OFF, modelPrint=OFF,
        contactPrint=OFF, historyPrint=OFF)
#endregion


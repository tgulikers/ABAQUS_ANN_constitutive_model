# database defining several material properties to use in abaqus

from abaqus import *


if allMats:
    for i in range(len(mdb.models.keys())):
        model_i = mdb.models[mdb.models.keys()[i]]

        model_i.Material(name='Steel1005')
        model_i.materials['Steel1005'].Elastic(table=((210000.0, 0.28),))
        model_i.materials['Steel1005'].Plastic(table=((350e9, 0.), (420e9, 0.15),))

        model_i.Material(name='alu6061')
        model_i.materials['alu6061'].Elastic(table=((68.9e3, 0.3),))
        model_i.materials['alu6061'].Plastic(table=((276.0, 0.), (310.0, 0.12),))

        model_i.Material(name='Composite_boyang1')
        model_i.materials['Composite_boyang1'].Elastic(table=((161000.0, 11000.0, 11000.0, 0.32, 0.32, 0.45,
                                                               5170.0, 5170.0, 3980.0),), type=ENGINEERING_CONSTANTS)
        model_i.materials['Composite_boyang1'].HashinDamageInitiation(
            table=((2800.0, 1700.0, 60.0, 125.0, 90.0, 90.0),))
        model_i.materials['Composite_boyang1'].hashinDamageInitiation.DamageEvolution(
            table=((100.0, 100.0, 0.22, 0.72),),
            type=ENERGY)
        model_i.materials['Composite_boyang1'].hashinDamageInitiation.DamageStabilization(fiberCompressiveCoeff=1e-05,
                                                                                          fiberTensileCoeff=1e-05,
                                                                                          matrixCompressiveCoeff=1e-05,
                                                                                          matrixTensileCoeff=1e-05)
else:
    currentModel.Material(name='SteelBoyang')
    currentModel.materials['SteelBoyang'].Elastic(table=((210000.0, 0.28),))
    currentModel.materials['SteelBoyang'].Plastic(table=((200.2, 0.), (246.0, 0.0235), (294.0, 0.0474), (374.0, 0.0935),
                                                        (437.0, 0.1377), (480.0, 0.1800)))

    currentModel.Material(name='alu6061')
    currentModel.materials['alu6061'].Elastic(table=((68.9e3, 0.3),))
    currentModel.materials['alu6061'].Plastic(table=((276.0, 0.), (310.0, 0.12),))

    currentModel.Material(name='Composite_boyang1')
    currentModel.materials['Composite_boyang1'].Elastic(table=((161000.0, 11000.0, 11000.0, 0.32, 0.32, 0.45,
                                                    5170.0, 5170.0, 3980.0),), type=ENGINEERING_CONSTANTS)
    currentModel.materials['Composite_boyang1'].HashinDamageInitiation(table=((2800.0, 1700.0, 60.0, 125.0, 90.0, 90.0),))
    currentModel.materials['Composite_boyang1'].hashinDamageInitiation.DamageEvolution(table=((100.0, 100.0, 0.22, 0.72),),
        type=ENERGY)
    currentModel.materials['Composite_boyang1'].hashinDamageInitiation.DamageStabilization(fiberCompressiveCoeff=1e-05,
        fiberTensileCoeff=1e-05, matrixCompressiveCoeff=1e-05, matrixTensileCoeff=1e-05)

    # week 4, course non-lin modelling!


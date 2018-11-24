The current folder contains the final Python code of the ANN module, developed for the master's thesis of:
T.H.E. Gulikers
Delft University of Technology
Date: 24-11-2018


The folder data_steelBrick_2D contains an example data set for a 2D steel plate, loaded in biaxial tension.
Run M_steelBrick_2D to train the ANN
The design of experiments is set up as follows for this data.
- Minimum load: strain = 0.0
- Maximum load: strain = 0.18
- Steps from minimum to maximum loads: 4
- Dimensions are fixed at height x width x depth = 10 x 10 x 1 mm


The folder data_compBrick_2D contains an example data set for a 2D composite plate with Hashin damage and an elliptical cut-out, loaded in biaxial tension and in-plane shear.
Run M_compBrick_2D to train the ANN
The design of experiments is set up as follows for this data.
- Minimum load: strain = 0.0
- Maximum load: strain = 0.02
- Steps from minimum to maximum loads: 4
- Dimensions are fixed at height x width x depth = 10 x 10 x 1 mm
- Lay-up: [45, 0, -45, 90]_s
- UD plies

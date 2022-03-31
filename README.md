# SGCo_IC
This work is associated with the paper "Assessing heterotopic searching strategy in hierarchical cosimulation for modeling the variables with inequality constraints" by Sultan Abulkhair and Nasser Madani published in Comptes Rendus. Géoscience. Repository belongs to Geostatistics research group in School of Mining and Geosciences, Nazarbayev University, Nur-Sultan, Kazakhstan. Please reference the original paper in case of using this code in your work:

Abulkhair, S., & Madani, N. (2021). Assessing heterotopic searching strategy in hierarchical cosimulation for modeling the variables with inequality constraints. Comptes Rendus. Géoscience, 353, 115-134. https://doi.org/10.5802/crgeos.58

# This file consists of 9 MATLAB functions:
1. sgco_ic.m - main code capable of implementing traditional sequential Gaussian, hierarchical and proposed cosimulation for variables 
with inequality constraint
2. create_paramfile_input.m - create default parameter file
3. cova.m - compute covariance values
4. setrot.m - set up matrix for rotation and reduction of coordinates
5. searchdata.m - search for neghboring data
6. weightcal.m - calculation of weights
7. backtr.m - back-transformation from Gaussian to original values
7. importfile.m - read data from GSLIB file
8. exportfile.m - print output into GSLIB file

Additionally, examples of target location Grid.out, files with variogram model vargfit (.cc, .nug and .mod), files with conversion tables (.trn) and parameter file for SGCo_IC are included.

# Instructions
1) Fill all input parameters in parameter file sgco_ic.par (to get default parameter file run create_paramfile_input.m).
2) Write following code in the MATLAB command line: sgco_ic('sgco_ic.par')
3) Output is GSLIB file that contains n relizations of primary and secondary variables.

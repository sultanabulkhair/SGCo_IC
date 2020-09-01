# SGCo_IC
This work is associated with the manuscript "Assessing Heterotopic Searching Strategy in Hierarchical Cosimulation Framework for Modeling the Variables with Inequality Constraints" by Sultan Abulkhair and Nasser Madani submitted to the Computers & Geosciences. Repository belongs to Geostatistics research group in School of Mining and Geosciences, Nazarbayev University, Nur-Sultan, Kazakhstan.
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

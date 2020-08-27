
function create_paramfile_input

%----------------------------------------------------
% Create a default parameter file for program SGCo-IC
%----------------------------------------------------
%

fid = fopen('sgco_ic.par','w');

fprintf(fid,'%1s\n','                  Parameters for SGCo_IC');
fprintf(fid,'%1s\n','                  **********************');
fprintf(fid,'%1s\n',' ');
fprintf(fid,'%1s\n','START OF PARAMETERS:');
fprintf(fid,'%1s\n','Grid.out                         % target location');
fprintf(fid,'%1s\n','1 2 3                            % columns for location coordinates');
fprintf(fid,'%1s\n','3 33 113 37                      % index of plane (0=all the planes), number of nodes ');
fprintf(fid,'%1s\n','nscore.out                       % file with conditioning data');
fprintf(fid,'%1s\n','1 2 3                            % columns for data coordinates');
fprintf(fid,'%1s\n','7 8                              % columns for Gaussian data values');
fprintf(fid,'%1s\n','Fe Al2O3                         % original variable names ');
fprintf(fid,'%1s\n','-10.0  10.0                      % trimming limits for Gaussian data');
fprintf(fid,'%1s\n','nscoreFe.trn                     % variable 1: file with conversion table (raw-Gaussian)');
fprintf(fid,'%1s\n','25.1 69.17                       %             minimum and maximum values for raw variable');
fprintf(fid,'%1s\n','1.0 1.0                          %             parameters for lower-tail and upper-tail extrapolation');
fprintf(fid,'%1s\n','nscoreAl2O3.trn                  % variable 2: file with conversion table (raw-Gaussian)');
fprintf(fid,'%1s\n','0.1  32.1                        %             minimum and maximum values for raw variable');
fprintf(fid,'%1s\n','1.0 1.0                          %             parameters for lower-tail and upper-tail extrapolation');
fprintf(fid,'%1s\n','0                                % Back-transformation? (0=yes, 1=no)');
fprintf(fid,'%1s\n','vargfit                          % basename for files with variogram models');
fprintf(fid,'%1s\n','2000 2000 2000                   % original data: maximum search radii in the rotated system');
fprintf(fid,'%1s\n','0 0 0                            % angles for search ellipsoid');
fprintf(fid,'%1s\n','0                                % divide into octants? 1=yes, 0=no');
fprintf(fid,'%1s\n','80                               % number of data per octant (if octant=1) or in total');
fprintf(fid,'%1s\n','2000 2000 2000                   % simulated nodes: maximum search radii in the rotated system');
fprintf(fid,'%1s\n','0 0 0                            % angles for search ellipsoid');
fprintf(fid,'%1s\n','0                                % divide into octants? 1=yes, 0=no');
fprintf(fid,'%1s\n','80                               % number of data per octant (if octant=1) or in total');
fprintf(fid,'%1s\n','1                                % co-simulation: 0=Conventional, 1=Hierarchical');
fprintf(fid,'%1s\n','0                                % if co-simulation=0 --> searching strategy: 0=single; 1=multiple');
fprintf(fid,'%1s\n','1                                % if co-simulation=1 --> 1th searching strategy: 0=single; 1=multiple');
fprintf(fid,'%1s\n','1                                % if co-simulation=1 --> 2th searching strategy: 0=single; 1=multiple');
fprintf(fid,'%1s\n','0                                % if co-simulation=1 --> inequality constraint: 0=yes; 1=no');
fprintf(fid,'%1s\n','-0.89 62                         % if inequality constraint=0: coefficients of linear restriction line: slope, intercept');
fprintf(fid,'%1s\n','100                              % number of realizations');
fprintf(fid,'%1s\n','9784498                          % seed for random number generation');
fprintf(fid,'%1s\n','sgco_ic.out                      % name of output file');
fprintf(fid,'%1s\n','3                                % number of decimals for values in the output file');
fprintf(fid,'%1s\n',' ');




fclose(fid);

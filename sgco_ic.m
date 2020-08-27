
function sgco_ic(simucoord,datacoord,datavalues,vname,limits,tableZY,zmin,zmax,tail,model,cc,nugget, ...
                 radius,angles,octant,ndata,ktype,nrealiz,nlines,ntok,seed,filename,nbdecimal);

%--------------------------------------------------------------------------------------------------
% Conditional co-simulation of cross-correlated Gaussian random fields via
% Sequential Gaussian co-Simulation
%--------------------------------------------------------------------------------------------------
%
% Author: Nasser Madani
% Last update: August 27, 2020
%
% INPUT:
%
%   CONDITIONING DATA
%   -----------------
%   datacoord    : data coordinates (n * 3 matrix; void for non-conditional simulations)
%   datavalues   : Gaussian conditioning data (n * nvar matrix; void for non-conditional simulations)
%   vname        : original variable names
%   limits       : trimming limits (inf and sup) for the Gaussian data (1 * 2 vector)
%
%   NORMAL SCORE TRANSFORMATION
%   ---------------------------
%   tableZY      : conversion tables between original and Gaussian values (void if no transformation is required)
%                   tableZY = [table_1;table_2;...;table_nvar]
%                   The first column of the table contains the original values, the second column their normal scores
%   zmin         : minimum values for the original variables (1 * nvar vector)
%   zmax         : maximum values for the original variables (1 * nvar vector)
%   tail         : parameters for lower-tail and upper-tail extrapolations (2 * nvar vector)
%
%   COREGIONALIZATION MODEL
%   -----------------------
%   model        : covariance model for the Gaussian random fields (nst * 7 matrix, where nst is the number of nested structures)
%                  Each row refers to a nested structure and is codified as: [type, ranges, angles]
%                  There are three ranges (along the rotated y, x and z axes) and three angles to define the coordinate rotation
%                  (azimuth, dip and plunge), see Deutsch and Journel, 1992, p. 25
%                  Available types:
%                    1: spherical
%                    2: exponential
%                    3: cubic
%                    4: Gaussian
%   cc           : sills or slopes of the nested structures (nst * nvar^2 matrix)
%   nugget       : nugget effect variance-covariance matrix (1 * nvar^2 vector)
%
%   SIMULATION PARAMETERS
%   ---------------------
%   ktype        : cokriging type: SCK=0, OCK=1
%   cross        : cross-validation? 1=yes, 0=no
%   nrealiz      : number of realizations to draw
%   nlines       : number of lines to use for each nested structure
%   ntok         : number of target points to simulate simultaneously
%   seed         : seed number for generating random values
%
%   OUTPUT OPTIONS
%   --------------
%   filename     : name of output file
%   nbdecimal    : number of decimals for the values in output file

%-----------------------------------------------------------------------------------------------------------------------------------

% It uses the following subroutines:
%   backtr.m                   : back-transformation from Gaussian to original values
%   cova.m                     : compute covariance values
%   create_paramfile_input.m  : create default parameter file
%   exportfile.m               : print output into GSLIB file
%   importfile.m               : read data from GSLIB file
%   setrot.m                   : set up matrix for rotation and reduction of coordinates
%   weightcal.m                : calculation of weights

%-----------------------------------------------------------------------------------------------------------

% User-defined parameters

warning('off','all');

%-----------------------------------------------------------------------------------------------------------------------------------

% Prompt for parameter file if no input is entered
%-------------------------------------------------

if (nargin < 1)
  disp('Which parameter file do you want to use?');
  paramfile = input('','s');
end


% If a single input is entered, it is the parameter file
%-------------------------------------------------------

if (nargin == 1)
  paramfile = simucoord;
end


% Read from parameter file
%-------------------------

if (nargin < 2)

  if isempty(paramfile)
    paramfile = 'sgco_ic.par';
  end

  fid = fopen(paramfile);

  if (fid < 0)
    disp('ERROR - The parameter file does not exist,');
    disp('        Check for the file and try again');
    disp(' ');
    disp('        creating a blank parameter file');
    disp(' ');
    disp('Stop - Program terminated.');
    create_paramfile_input;
    return;
  else
    disp(' ');
    disp('Reading parameter file');
  end

  % The parameter file does exist
  fgets(fid); fgets(fid); fgets(fid); fgets(fid);


  tline0 = fgets(fid);
  i = (find(tline0 == ' '));
  if ~isempty(i), tline0 = tline0(1:i(1)-1); end

  tline = fgets(fid);
  index = str2num(tline);

  [filetitle1,nfield1,header1,grid] = importfile(tline0,[],1,index);
  
  tline = fgets(fid);
  plane = str2num(tline);
 

   tline0 = fgets(fid);
  i = (find(tline0 == ' '));
  if ~isempty(i), tline0 = tline0(1:i(1)-1); end

  tline = fgets(fid);
  ind_coord_data = str2num(tline);

  tline = fgets(fid);
  ind_ps = str2num(tline);

  [filetitle2,nfield2,header2,datacoord,data] = importfile(tline0,[],1,ind_coord_data,ind_ps);
 
    
  nvar=2;
  
  tline = fgets(fid);
  vname = tline;

  tline = fgets(fid);
  limits = str2num(tline);

  
  tableZY = zeros(0,2);
  zmin = -inf*ones(1,nvar);
  zmax = inf*ones(1,nvar);
  tail = ones(2,nvar);
  for j = 1:nvar
    tline = fgets(fid);
    i = (find(tline == ' '));
    if ~isempty(i), tline = tline(1:i(1)-1); end
    fid2(j) = fopen(tline);
    if (fid2(j) > -1)
      tablej = dlmread(tline);
      lj = size(tablej,1);
      J = find( (tablej(2:lj-1,1)==tablej(1:lj-2,1)) & (tablej(2:lj-1,1)==tablej(3:lj,1)) );
      tablej(J+1,:) = [];
      tableZY = [tableZY;tablej(:,1:2)];
    end
    tline = fgets(fid);
    tline = str2num(tline);
    zmin(j) = tline(1); zmax(j) = tline(2);
    tline = fgets(fid);
    tail(:,j) = str2num(tline)';
  end

  tline = fgets(fid);
  back = str2num(tline);
  
  
  
  tline = fgets(fid);
  j = find((tline == ' ')|(tline == char(10))|(tline == char(13)));
  if ~isempty(j), tline = tline(1:j(1)-1); end
  fid3 = fopen([tline,'.mod']);
  if (fid3 > -1), model = load([tline,'.mod'],'-ascii'); fclose(fid3); else, error(['Unable to read file ',tline,'.mod']); end
  fid3 = fopen([tline,'.cc']);
  if (fid3 > -1), cc = load([tline,'.cc'],'-ascii'); fclose(fid3); else, error(['Unable to read file ',tline,'.cc']); end
  fid3 = fopen([tline,'.nug']);
  if (fid3 > -1), nugget = load([tline,'.nug'],'-ascii'); fclose(fid3); else, error(['Unable to read file ',tline,'.nug']); end

  tline = fgets(fid);
  search_radii_data = str2num(tline);

  tline = fgets(fid);
  angle_data = str2num(tline);
  
  tline = fgets(fid);
  divide_data = str2num(tline);

  tline = fgets(fid);
  noctant_data = str2num(tline);

  
  tline = fgets(fid);
  search_radii_simul = str2num(tline);

  tline = fgets(fid);
  angle_simul = str2num(tline);
  
  tline = fgets(fid);
  divide_simul = str2num(tline);

  tline = fgets(fid);
  noctant_simul = str2num(tline);
  
  tline = fgets(fid);
  cosim = str2num(tline);
  
  tline = fgets(fid);
  cosim_1 = str2num(tline);
  
  tline = fgets(fid);
  cosim_2 = str2num(tline);
  
  tline = fgets(fid);
  cosim_3 = str2num(tline);
  
  tline = fgets(fid);
  ic = str2num(tline);
  
  tline = fgets(fid);
  coefficient = str2num(tline);
  
  tline = fgets(fid);
  nrealiz = str2num(tline);
  
  tline = fgets(fid);
  seed = str2num(tline);

  filename = fgets(fid);
  i = (find(filename == ' '));
  if ~isempty(i), filename = filename(1:i(1)-1); end

  tline = fgets(fid);
  nbdecimal = str2num(tline);

  fclose(fid);
  for j = 1:1, if fid2(j)>-1, fclose(fid2(j)); end; end

end


disp(' ');
disp('Preparing simulation');

% Define default values
%----------------------

ind_coord = [1:3]; % index of coordinates [x Y Z]
ngrid = size(grid,1);

nx = plane(2); ny=plane(3);nz=plane(4);
if plane(1)>0
    index = 1:ngrid;
    index = reshape(index,nx*ny,nz);
    grid = grid(index(:,plane(1)),:);
end

ndata=size(data,1);
ngrid = size(grid,1);
coord_g = grid(:,1:3);
coord_d = datacoord;



data_1(:,1) = data(:,1);
data_1(:,2) = data(:,2);

search_rotationmatrix_data = setrot([1 search_radii_data angle_data],1);
search_rotationmatrix_simul = setrot([1 search_radii_simul angle_simul],1);
simul_temp = ones(ngrid,2*nrealiz).*-99;
simul = simul_temp;
index_seq = randsample(ngrid,ngrid);



nst = size(model,1);     % number of nested structures
nvar = sqrt(size(cc,2)); % number of variables
progress = 0;
for j = 1:nvar
  while vname(1)==' ', vname(1) = []; end
  i = find(vname==' ');
  names{j} = vname(1:i(1)-1);
  vname(1:i(1)-1) = [];
end

if nvar > floor(nvar), error('The number of columns in the sill matrix (cc) is inconsistent'); end

sill = zeros(nvar,nvar,nst);
for i = 1:nst
  sill(:,:,i) = reshape(cc(i,:),nvar,nvar);
  R = setrot(model,i);
  model_rotationmatrix(:,:,i) = R;
  if max(abs(sill(:,:,i)-sill(:,:,i)'))>100*eps, error(['The sill matrix for structure n?',num2str(i),' is not symmetric']); end
  [eigenvectors,eigenvalues] = eig(sill(:,:,i));
  if min(diag(eigenvalues))<0, error(['The sill matrix for structure n?',num2str(i),' is not positive semi-definite']); end
end

sillnugget = reshape(nugget,nvar,nvar);
if max(abs(sillnugget-sillnugget'))>100*eps, error(['The sill matrix for the nugget effect is not symmetric']); end
[eigenvectors,eigenvalues] = eig(sillnugget);
A0 = sqrt(eigenvalues)*eigenvectors';
if min(diag(eigenvalues))<0, error(['The sill matrix for the nugget effect is not positive semi-definite']); end


if ~isempty(tableZY)
  tableZY = tableZY(:,1:2);
  length_table = size(tableZY,1);
  index_table = [0;find(tableZY(2:length_table,2)<tableZY(1:length_table-1,2));length_table];
  tableZY(:,2) = tableZY(:,2) + 1e-10*[1:length_table]'; % avoid tied values in the conversion table
end


% Create seed numbers
%--------------------

rand('state',seed);
randn('state',seed);



if cosim <1  % Conventional co-simulation
for j=1:ngrid
    i = index_seq(j);
    simul_temp = simul;
   
    
    if cosim_1==0 % Single search-----
        [I_d] = searchdata(coord_d,data_1,coord_g(i,1:3),search_rotationmatrix_data,divide_data,noctant_data,1e-20);
        I_d_1=I_d;
        I_d_2=I_d;
        index_mising_1 = find(data_1(I_d_1,1)==-99);
        I_d_1(index_mising_1)=[];
        datacoord_d_1 = coord_d(I_d_1,:);
        datavalues_d_1 = data_1(I_d_1,1);
        index_mising_2 = find(data_1(I_d_2,2)==-99);
        I_d_2(index_mising_2)=[];
        datacoord_d_2 = coord_d(I_d_2,:);
        datavalues_d_2 = data_1(I_d_2,2);
        
        [I_g] = searchdata(coord_g,simul_temp,coord_g(i,1:3),search_rotationmatrix_simul,divide_simul,noctant_simul,1e-20);
        I_g_1=I_g;
        I_g_2=I_g;
        index_mising_3 = find(simul_temp(I_g_1,1)==-99);
        I_g_1(index_mising_3)=[];
        datacoord_g_1 = coord_g(I_g_1,:);
        datavalues_g_1 = simul_temp(I_g_1,1:2:2*nrealiz); 
    
        index_mising_4 = find(simul_temp(I_g_2,2)==-99);
        I_g_2(index_mising_4)=[];
        datacoord_g_2 = coord_g(I_g_2,:);
        datavalues_g_2 = simul_temp(I_g_2,2:2:2*nrealiz); 
          
        index_v1 (1) = size(I_d_1,1);  index_v1 (2) = size(I_g_1,1);
        index_v2 (1) = size(I_d_2,1);  index_v2 (2) = size(I_g_2,1);
        
        coord_1 = [datacoord_d_1; datacoord_g_1; datacoord_d_2; datacoord_g_2;  coord_g(i,1:3)];
        [lambda,sigma] = weightcal(coord_1,model,cc,nugget,ktype,index_v1,index_v2);
        data_2 = [kron(datavalues_d_1(:,1),ones(1,nrealiz));datavalues_g_1;kron(datavalues_d_2(:,1),ones(1,nrealiz));datavalues_g_2];
   
    elseif cosim_1==1 % Multiple search-----
        ktype=0;
        % Searching the first variable on conditioning data
        index_mising_1 = find(data_1(:,1)==-99);
        data_1_1= data_1(:,1);
        coord_d_1=coord_d;
        data_1_1 (index_mising_1)=[];
        coord_d_1(index_mising_1,:)=[];
        
        [I_d_1] = searchdata(coord_d_1,data_1_1,coord_g(i,1:3),search_rotationmatrix_data,divide_data,noctant_data,1e-20);
        datacoord_d_1 = coord_d_1(I_d_1,:);
        datavalues_d_1 = data_1_1(I_d_1,:);
        
         % Searching the second variable on conditioning data
        index_mising_2 = find(data_1(:,2)==-99);
        data_1_2= data_1(:,2);
        coord_d_2=coord_d;
        data_1_2 (index_mising_2)=[];
        coord_d_2(index_mising_2,:)=[];
          
        [I_d_2] = searchdata(coord_d_2,data_1_2,coord_g(i,1:3),search_rotationmatrix_data,divide_data,noctant_data,1e-20);
        datacoord_d_2 = coord_d_2(I_d_2,:);
        datavalues_d_2 = data_1_2(I_d_2,:);
         
        % Searching the first variable on simulated values
        index_mising_3 = find(simul_temp(:,1)==-99);
        simul_temp_1 = simul_temp(:,1:2:2*nrealiz);
        coord_g_1 = coord_g;
        simul_temp_1(index_mising_3,:)=[];
        coord_g_1(index_mising_3,:)=[];
          
        [I_g_1] = searchdata(coord_g_1,simul_temp_1,coord_g(i,1:3),search_rotationmatrix_simul,divide_simul,noctant_simul,1e-20);
        datacoord_g_1 = coord_g_1(I_g_1,:);
        datavalues_g_1 = simul_temp_1(I_g_1,:);  
    
        % Searching the second variable on simulated values
        index_mising_4 = find(simul_temp(:,2)==-99);
        simul_temp_2 = simul_temp(:,2:2:2*nrealiz);
        coord_g_2 = coord_g;
        simul_temp_2(index_mising_4,:)=[];
        coord_g_2(index_mising_4,:)=[];
         
        [I_g_2] = searchdata(coord_g_2,simul_temp_2,coord_g(i,1:3),search_rotationmatrix_simul,divide_simul,noctant_simul,1e-20);
        datacoord_g_2 = coord_g_2(I_g_2,:);
        datavalues_g_2 = simul_temp_2(I_g_2,:);  
        index_v1 (1) = size(I_d_1,1);  index_v1 (2) = size(I_g_1,1);
        index_v2 (1) = size(I_d_2,1);  index_v2 (2) = size(I_g_2,1);
        
         coord_1 = [datacoord_d_1; datacoord_g_1; datacoord_d_2; datacoord_g_2;  coord_g(i,1:3)];
        [lambda,sigma] = weightcal(coord_1,model,cc,nugget,ktype,index_v1,index_v2);
        data_2 = [kron(datavalues_d_1,ones(1,nrealiz));datavalues_g_1;kron(datavalues_d_2,ones(1,nrealiz));datavalues_g_2];
        
       
    end
    
    [V,D] = eig(sigma);
    sqrtsigma = V*sqrt(max(D,0));
     estimation = (data_2'*lambda)';
    simul(i,:) = reshape(estimation + sqrtsigma*randn(2,nrealiz),1,2*nrealiz);
      
     % Report on progress from time to time
    progress2 = 10*floor(10*j/ngrid);
    if (progress2 > progress)
      disp(['  ',num2str(progress2),'% completed']);
      progress = progress2;
      pause(0.001);
    end
    
end




elseif cosim<2 % Hierarchical co-simulation
   % First run for simulating the first variable by FC
    disp(' ');
    disp('First run for simulating the first variable');
 for j=1:ngrid
    i = index_seq(j);
    simul_temp = simul;
   
    
    if cosim_2==0 % Single search-----
        ktype=0;
        [I_d] = searchdata(coord_d,data_1,coord_g(i,1:3),search_rotationmatrix_data,divide_data,noctant_data,1e-20);
        I_d_1=I_d;
        I_d_2=I_d;
        index_mising_1 = find(data_1(I_d_1,1)==-99);
        I_d_1(index_mising_1)=[];
        datacoord_d_1 = coord_d(I_d_1,:);
        datavalues_d_1 = data_1(I_d_1,1);
     
        index_mising_2 = find(data_1(I_d_2,2)==-99);
        I_d_2(index_mising_2)=[];
        datacoord_d_2 = coord_d(I_d_2,:);
        datavalues_d_2 = data_1(I_d_2,2);
          
        [I_g] = searchdata(coord_g,simul_temp,coord_g(i,1:3),search_rotationmatrix_simul,divide_simul,noctant_simul,1e-20);
        I_g_1=I_g;
        I_g_2=I_g;
        index_mising_3 = find(simul_temp(I_g_1,1)==-99);
        I_g_1(index_mising_3)=[];
        datacoord_g_1 = coord_g(I_g_1,:);
        datavalues_g_1 = simul_temp(I_g_1,1:2:2*nrealiz); 
    
        index_mising_4 = find(simul_temp(I_g_2,2)==-99);
        I_g_2(index_mising_4)=[];
        datacoord_g_2 = coord_g(I_g_2,:);
        datavalues_g_2 = simul_temp(I_g_2,2:2:2*nrealiz); 
           
        index_v1 (1) = size(I_d_1,1);  index_v1 (2) = size(I_g_1,1);
        index_v2 (1) = size(I_d_2,1);  index_v2 (2) = size(I_g_2,1);
        
        coord_1 = [datacoord_d_1; datacoord_g_1; datacoord_d_2; datacoord_g_2;  coord_g(i,1:3)];
        [lambda,sigma] = weightcal(coord_1,model,cc,nugget,ktype,index_v1,index_v2);
        data_2 = [kron(datavalues_d_1(:,1),ones(1,nrealiz));datavalues_g_1;kron(datavalues_d_2(:,1),ones(1,nrealiz));datavalues_g_2];
         
    elseif cosim_2==1 % Multiple search-----
        ktype=0;
        % Searching the first variable on conditioning data
        index_mising_1 = find(data_1(:,1)==-99);
        data_1_1= data_1(:,1);
        coord_d_1=coord_d;
        data_1_1 (index_mising_1)=[];
        coord_d_1(index_mising_1,:)=[];
        
        [I_d_1] = searchdata(coord_d_1,data_1_1,coord_g(i,1:3),search_rotationmatrix_data,divide_data,noctant_data,1e-20);
        datacoord_d_1 = coord_d_1(I_d_1,:);
        datavalues_d_1 = data_1_1(I_d_1,:);
        
         % Searching the second variable on conditioning data
        index_mising_2 = find(data_1(:,2)==-99);
        data_1_2= data_1(:,2);
        coord_d_2=coord_d;
        data_1_2 (index_mising_2)=[];
        coord_d_2(index_mising_2,:)=[];
          
        [I_d_2] = searchdata(coord_d_2,data_1_2,coord_g(i,1:3),search_rotationmatrix_data,divide_data,noctant_data,1e-20);
        datacoord_d_2 = coord_d_2(I_d_2,:);
        datavalues_d_2 = data_1_2(I_d_2,:);
         
        % Searching the first variable on simulated values
        index_mising_3 = find(simul_temp(:,1)==-99);
        simul_temp_1 = simul_temp(:,1:2:2*nrealiz);
        coord_g_1 = coord_g;
        simul_temp_1(index_mising_3,:)=[];
        coord_g_1(index_mising_3,:)=[];
     
        [I_g_1] = searchdata(coord_g_1,simul_temp_1,coord_g(i,1:3),search_rotationmatrix_simul,divide_simul,noctant_simul,1e-20);
        datacoord_g_1 = coord_g_1(I_g_1,:);
        datavalues_g_1 = simul_temp_1(I_g_1,:);  
    
        % Searching the second variable on simulated values
      
        index_mising_4 = find(simul_temp(:,2)==-99);
        simul_temp_2 = simul_temp(:,2:2:2*nrealiz);
        coord_g_2 = coord_g;
        simul_temp_2(index_mising_4,:)=[];
        coord_g_2(index_mising_4,:)=[];
          
        [I_g_2] = searchdata(coord_g_2,simul_temp_2,coord_g(i,1:3),search_rotationmatrix_simul,divide_simul,noctant_simul,1e-20);
        datacoord_g_2 = coord_g_2(I_g_2,:);
        datavalues_g_2 = simul_temp_2(I_g_2,:);  
         
         index_v1 (1) = size(I_d_1,1);  index_v1 (2) = size(I_g_1,1);
         index_v2 (1) = size(I_d_2,1);  index_v2 (2) = size(I_g_2,1);
        
        coord_1 = [datacoord_d_1; datacoord_g_1; datacoord_d_2; datacoord_g_2;  coord_g(i,1:3)];
        [lambda,sigma] = weightcal(coord_1,model,cc,nugget,ktype,index_v1,index_v2);
        data_2 = [kron(datavalues_d_1,ones(1,nrealiz));datavalues_g_1;kron(datavalues_d_2,ones(1,nrealiz));datavalues_g_2];
    
    end
    
   
    [V,D] = eig(sigma);
    sqrtsigma = V*sqrt(max(D,0));
    estimation = (data_2'*lambda)';
    simul(i,:) = reshape(estimation + sqrtsigma*randn(2,nrealiz),1,2*nrealiz);

    
    
       % Report on progress from time to time
    progress2 = 10*floor(10*j/ngrid);
    if (progress2 > progress)
      disp(['  ',num2str(progress2),'% completed']);
      progress = progress2;
      pause(0.001);
    end
    
 end
    
  

 
 
 
 simul (:,2:2:2*nrealiz)=-99;

 a=zeros(size(simul,1),2*nrealiz);
 table_2 = tableZY(index_table(2)+1:index_table(2+1),:);
 progress = 0;

 for m=1:1 
   % Second run for simulating the second variable by Different searchnig
   % strategy
       disp(' ');
    disp('Second run for simulating the second variable');
 for j=1:ngrid
    i = index_seq(j);
    simul_temp = simul;
   
    
      if cosim_3==0 % Single search cokriging
            ktype = 5;
           % Searching the first variable on conditioning data
        [I_d] = searchdata(coord_d,data_1,coord_g(i,1:3),search_rotationmatrix_data,divide_data,noctant_data,1e-20);
        I_d_1=I_d;
        I_d_2=I_d;
        index_mising_1 = find(data_1(I_d_1,1)==-99);
        I_d_1(index_mising_1)=[];
        datacoord_d_1 = coord_d(I_d_1,:);
        datavalues_d_1 = data_1(I_d_1,1);
              
        index_mising_2 = find(data_1(I_d_2,2)==-99);
        I_d_2(index_mising_2)=[];
        datacoord_d_2 = coord_d(I_d_2,:);
        datavalues_d_2 = data_1(I_d_2,2);
         
        [I_g] = searchdata(coord_g,simul_temp,coord_g(i,1:3),search_rotationmatrix_simul,divide_simul,noctant_simul,1e-20);
        I_g_1=I_g;
        I_g_2=I_g;
        index_mising_3 = find(simul_temp(I_g_1,1)==-99);
        I_g_1(index_mising_3)=[];
        datacoord_g_1 = coord_g(I_g_1,:);
        datavalues_g_1 = simul_temp(I_g_1,1:2:2*nrealiz); 
    
        index_mising_4 = find(simul_temp(I_g_2,2)==-99);
        I_g_2(index_mising_4)=[];
        datacoord_g_2 = coord_g(I_g_2,:);
        datavalues_g_2 = simul_temp(I_g_2,2:2:2*nrealiz); 
    
        index_v1 (1) = size(I_d_1,1);index_v1 (2) = size(I_g_1,1);
        index_v2 (1) = size(I_d_2,1);  index_v2 (2) = size(I_g_2,1);
        
        coord_1 = [datacoord_d_2; datacoord_g_2;datacoord_d_1; datacoord_g_1;  coord_g(i,1:3)]; 
        [lambda,sigma] = weightcal(coord_1,model,cc,nugget,ktype,index_v1,index_v2);
        data_2 = [kron(datavalues_d_2,ones(1,nrealiz));datavalues_g_2;...
        kron(datavalues_d_1,ones(1,nrealiz));datavalues_g_1;simul(i,1:2:2*nrealiz)];

        estimation = (data_2'*lambda)';
  
        if ic==0
              a=estimation(1,:) + sqrt(max(0,sigma(1)))*randn(1,nrealiz);
              simul_final = a;
              simul_b=simul(i,1:2:2*nrealiz);
              if ~isempty(tableZY)
                     simul_b(:,1:nrealiz) = backtr(simul(i,1:2:2*nrealiz),tableZY(index_table(1)+1:index_table(1+1),:),zmin(1),zmax(1),tail(:,1));
              end
               if ~isempty(tableZY)
                     simul_b_2(:,1:nrealiz) = backtr(a,tableZY(index_table(2)+1:index_table(2+1),:),zmin(2),zmax(2),tail(:,2));

              end
               maxu = coefficient(1,1)*simul_b+coefficient(1,2);
               maxu_n=backtr(maxu,[table_2(:,2) table_2(:,1)],min(table_2(:,2)),max(table_2(:,2)),[1 1]);
              
               index = a<=maxu_n;
               
               index_0 = find(index==0);
               
               if ~isempty(index_0)
                   maxu_n=backtr(maxu(:,index_0),[table_2(:,2) table_2(:,1)],min(table_2(:,2)),max(table_2(:,2)),[1 1]);
                   minu_n = -inf;
                   minu_n=(minu_n-estimation(1,index_0))./sqrt(max(0,sigma(1)));
                   maxu_n=(maxu_n-estimation(1,index_0))./sqrt(max(0,sigma(1)));
                   minu_n = normcdf(minu_n);
                   maxu_n = normcdf(maxu_n);
                   v = norminv(minu_n + (maxu_n-minu_n).*rand(1,size(index_0,2)));
                   b =estimation(1,index_0) + sqrt(max(0,sigma(1))).*v;              
                   simul_final(:,index_0)=b;
                   
               end
              simul(i,2:2:2*nrealiz)=simul_final;
                 else 
              simul(i,2:2:2*nrealiz) = estimation(1,:) + sqrt(max(0,sigma(1)))*randn(1,nrealiz);
           end 
      
       
          elseif cosim_3==1 % Multiple search cokriging  
            ktype = 5;
       % Searching the first variable on conditioning data
        index_mising_1 = find(data_1(:,1)==-99);
        data_1_1= data_1(:,1);
        coord_d_1=coord_d;
        data_1_1 (index_mising_1)=[];
        coord_d_1(index_mising_1,:)=[];
        
        [I_d_1] = searchdata(coord_d_1,data_1_1,coord_g(i,1:3),search_rotationmatrix_data,divide_data,noctant_data,1e-20);
        datacoord_d_1 = coord_d_1(I_d_1,:);
        datavalues_d_1 = data_1_1(I_d_1,:);
        
         % Searching the second variable on conditioning data
        index_mising_2 = find(data_1(:,2)==-99);
        data_1_2= data_1(:,2);
        coord_d_2=coord_d;
        data_1_2 (index_mising_2)=[];
        coord_d_2(index_mising_2,:)=[];
          
        [I_d_2] = searchdata(coord_d_2,data_1_2,coord_g(i,1:3),search_rotationmatrix_data,divide_data,noctant_data,1e-20);
        datacoord_d_2 = coord_d_2(I_d_2,:);
        datavalues_d_2 = data_1_2(I_d_2,:);
          
         % Searching the first variable on simulated values
        index_mising_3 = find(simul_temp(:,1)==-99);
        simul_temp_1 = simul_temp(:,1:2:2*nrealiz);
        coord_g_1 = coord_g;
        simul_temp_1(index_mising_3,:)=[];
        coord_g_1(index_mising_3,:)=[];
          
        [I_g_1] = searchdata(coord_g_1,simul_temp_1,coord_g(i,1:3),search_rotationmatrix_simul,divide_simul,noctant_simul,1e-20);
        datacoord_g_1 = coord_g_1(I_g_1,:);
        datavalues_g_1 = simul_temp_1(I_g_1,:);  
    
        % Searching the second variable on simulated values
        index_mising_4 = find(simul_temp(:,2)==-99);
        simul_temp_2 = simul_temp(:,2:2:2*nrealiz);
        coord_g_2 = coord_g;
        simul_temp_2(index_mising_4,:)=[];
        coord_g_2(index_mising_4,:)=[];
           
        [I_g_2] = searchdata(coord_g_2,simul_temp_2,coord_g(i,1:3),search_rotationmatrix_simul,divide_simul,noctant_simul,1e-20);
        datacoord_g_2 = coord_g_2(I_g_2,:);
        datavalues_g_2 = simul_temp_2(I_g_2,:);  
      
        index_v1 (1) = size(I_d_1,1);  index_v1 (2) = size(I_g_1,1);
        index_v2 (1) = size(I_d_2,1);  index_v2 (2) = size(I_g_2,1);
            
        coord_1 = [datacoord_d_2; datacoord_g_2; datacoord_d_1; datacoord_g_1; coord_g(i,1:3)]; 
        [lambda,sigma] = weightcal(coord_1,model,cc,nugget,ktype,index_v1,index_v2);
        data_2 = [kron(datavalues_d_2,ones(1,nrealiz));datavalues_g_2;...
        kron(datavalues_d_1,ones(1,nrealiz));datavalues_g_1;simul(i,1:2:2*nrealiz)];
        estimation = (data_2'*lambda)';
       
        if ic==0
              a=estimation(1,:) + sqrt(max(0,sigma(1)))*randn(1,nrealiz);
              simul_final = a;
              simul_b=simul(i,1:2:2*nrealiz);
              if ~isempty(tableZY)
                     simul_b(:,1:nrealiz) = backtr(simul(i,1:2:2*nrealiz),tableZY(index_table(1)+1:index_table(1+1),:),zmin(1),zmax(1),tail(:,1));

              end
               if ~isempty(tableZY)
                     simul_b_2(:,1:nrealiz) = backtr(a,tableZY(index_table(2)+1:index_table(2+1),:),zmin(2),zmax(2),tail(:,2));
              end
               maxu = coefficient(1,1)*simul_b+coefficient(1,2);
               maxu_n=backtr(maxu,[table_2(:,2) table_2(:,1)],min(table_2(:,2)),max(table_2(:,2)),[1 1]);
              
               index = a<=maxu_n;
               
               index_0 = find(index==0);
               
               if ~isempty(index_0)
                   maxu_n=backtr(maxu(:,index_0),[table_2(:,2) table_2(:,1)],min(table_2(:,2)),max(table_2(:,2)),[1 1]);
                   minu_n =-inf;
                   minu_n=(minu_n-estimation(1,index_0))./sqrt(max(0,sigma(1)));
                   maxu_n=(maxu_n-estimation(1,index_0))./sqrt(max(0,sigma(1)));
                   minu_n = normcdf(minu_n);
                   maxu_n = normcdf(maxu_n);
                   v = norminv(minu_n + (maxu_n-minu_n).*rand(1,size(index_0,2)));
                   b =estimation(1,index_0) + sqrt(max(0,sigma(1))).*v;              
                   simul_final(:,index_0)=b;
                   
               end
                simul(i,2:2:2*nrealiz)=simul_final;
         else 
              simul(i,2:2:2*nrealiz) = estimation(1,:) + sqrt(max(0,sigma(1)))*randn(1,nrealiz);
        
        end

              
      end
       
    
      % Report on progress from time to time
    progress2 = 10*floor(10*j/ngrid);
    if (progress2 > progress)
      disp(['  ',num2str(progress2),'% completed']);
      progress = progress2;
      pause(0.001);
    end
    
 end   
 
end
end



 % Back transform to original scale
      %---------------------------------
if back==0
    if ~isempty(tableZY)
        for i = 1:nvar
        simul(:,i:2:2*nrealiz) = backtr(simul(:,i:2:2*nrealiz),tableZY(index_table(i)+1:index_table(i+1),:),zmin(i),zmax(i),tail(:,i));
        end
    end
end
 
% Write in output file
%---------------------

filetitle = 'Simulated values';
nfield = nrealiz*nvar;
k = 1;
for i = 1:nrealiz
  for j = 1:nvar
    header{k} = ['Realization no.',int2str(i),' - ',names{j}];
    k = k+1;
  end
end
exportfile(filename,filetitle,nfield,header,reshape(simul,ngrid,nvar*nrealiz),[],nbdecimal(1),1,1);




  % CLOSE THE OUTPUT FILE
  %----------------------

 % fclose(fid); 
    
      
disp(' ')
disp('... Done...')
end

  




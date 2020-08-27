function [filetitle,nfield,header,varargout] = importfile(filename,format,numeric,varargin);

%------------------------------
% Import file with GSLIB header
%------------------------------
%
%
% INPUT:
%   filename : name for input file
%   format   : format for file columns ('%f' = float; '%d'=integer; '%s'=string)
%   numeric  : import as numeric data? (1=yes, 0=string data)
%   varargin : index of columns to extract

%-----------------------------------------------------------------------------------------------

warning('off','all');

fid = fopen(filename);
if (fid < 0)
  filetitle = [];
  nfield = 0;
  header = [];
  for i = 1:nargout-3
    varargout{i} = [];
  end
  disp(['  Warning: cannot import ',filename]);
  return
end

tline = fgetl(fid);
while tline(1)==[' '], tline(1) = []; end
while tline(length(tline))==[' '], tline(length(tline)) = []; end
filetitle = tline;

tline = fgetl(fid);
nfield = str2num(tline);

for i = 1:nfield(1)
  tline = fgetl(fid);
  while tline(1)==[' '], tline(1) = []; end
  while tline(length(tline))==[' '], tline(length(tline)) = []; end
  header{i} = tline;
end

if nargin<2, format=[]; end
if isempty(format)
  for j = 1:nfield(1), format = [format,'%f ']; end
end
if isnan(format)
  format = [];
  for j = 1:nfield(1), format = [format,'%s ']; end
end
allinput = textscan(fid,format);

if nargin<3, numeric = 1; end
if isempty(numeric), numeric = 1; end

if (nargin<4)&&(numeric==0)
  varargout{1} = allinput;
  return;
end

if (nargin<4)
  T = [];
  for j = 1:nfield(1)
    T = [T allinput{j}];
  end
  varargout{1} = T;
  return;
end

for i = 1:nargin-3
  index = varargin{i};
  if isempty(index)
    varargout{i} = [];
    continue;
  end
  index0 = find(index==0);
  if length(index0) == length(index)
    varargout{i} = [];
    continue;
  end
  T = [];
  for j = 1:length(index)
    if index(j)==0
      T = [T ones(size(allinput{1},1),1)];
    else
      T = [T allinput{index(j)}];
    end
  end
  varargout{i} = T;
end

fclose(fid);

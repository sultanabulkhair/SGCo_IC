function fid = exportfile(filename,filetitle,nfield,header,filecontent,format,nbdecimal,fileopen,fileclose);

%------------------------------
% Export file with GSLIB header
%------------------------------
%
%
% INPUT:
%   filename    : name for output file, or file identifier if fileopen=0
%   filetitle   : header title
%   nfield      : number of fields
%   header      : header
%   filecontent : contents to export
%   format      : format for file columns ('%f' = float; '%d'=integer; '%s'=string)
%   nbdecimal   : number of decimals for float outputs
%   fileopen    : open output file? 1=yes, 0=no
%   fileclose   : close output file? 1=yes, 0=no

%-----------------------------------------------------------------------------------------------

warning('off','all');

if nargin<9, fileclose = 1; end
if nargin<8, fileopen = 1; end
if nargin<7, nbdecimal = 5; end
if nargin<6, format=[]; end

if fileopen > 0
  fid = fopen(filename,'w');
  fprintf(fid,'%1s\n',filetitle);
  fprintf(fid,'%1s\n',int2str(nfield));
  for i = 1:nfield(1)
    fprintf(fid,'%1s\n',header{i});
  end
else
  fid = filename;
end

if isempty(format)
  for i = 1:nfield(1), format = [format,'  %.',int2str(nbdecimal(1)),'f']; end
end

if iscell(filecontent)
  for i = 1:size(filecontent{1},1)
    clear A;
    for j = 1:size(filecontent,2)
      if iscell(filecontent{j}(i))
        tempcell = filecontent{j}(i);
        A{1,j} = tempcell{1};
      else
        A{1,j} = filecontent{j}(i);
      end
    end
    fprintf(fid,[format,'\n'],A{1,:});
  end
elseif ~isempty(filecontent)
  I = find(isnan(filecontent(:)));
  filecontent(I) = -99*ones(size(I));
  fprintf(fid,[format,'\n'],filecontent');
end

if fileclose > 0
  fclose(fid);
end

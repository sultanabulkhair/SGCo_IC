
function z = backtr(y,table,zmin,zmax,tail)

%----------------------------------------------------
% Back-transformation from Gaussian to original scale
%----------------------------------------------------
%

p = size(table,1);
z = y;
if isempty(table), return; end


% Values in the lower tail (exponential extrapolation)
%-----------------------------------------------------

z1 = table(1,1);
y1 = table(1,2);
I1 = find(y(:)<y1);

if ~isempty(I1)
  b0 = (z1-zmin)*exp(-tail(1)*y1);
  z(I1) = zmin + b0*exp(tail(1)*y(I1));
end


% Values in the upper tail (exponential extrapolation)
%-----------------------------------------------------

zp = table(p,1);
yp = table(p,2);
I2 = find(y(:)>yp);

if ~isempty(I2)
  bp = (zp-zmax)*exp(tail(2)*yp);
  z(I2) = zmax + bp*exp(-tail(2)*y(I2));
end


% Within-class values
%--------------------

I = [I1(:);I2(:)];
I3 = [1:length(y(:))]';
I3(I) = [];

if ~isempty(I3)
  z(I3) = interp1(table(:,2)+1e-12*[1:size(table,1)]',table(:,1),y(I3)); % Table lookup
end

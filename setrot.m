
function rotred_matrix = setrot(model,it);

%------------------------------------------------------------
% Set up the matrix to transform the Cartesian coordinates to
% coordinates that account for angles and anisotropy
%------------------------------------------------------------
%
%
% INPUT:
%   model: nested variogram model
%   it: index of the nested structure
%   The rotation is performed according to the GSLIB conventions

deg2rad = pi/180;
ranges = model(it,2:4);
angles = model(it,5:7);

% matrix of coordinate reduction
redmat = diag(1./(eps+ranges));

a = (90-angles(1))*deg2rad;
b = -angles(2)*deg2rad;
c = angles(3)*deg2rad;

cosa = cos(a);
sina = sin(a);
cosb = cos(b);
sinb = sin(b);
cosc = cos(c);
sinc = sin(c);

rotmat = zeros(3,3);
rotmat(1,1) = cosb * cosa;
rotmat(1,2) = cosb * sina;
rotmat(1,3) = -sinb;
rotmat(2,1) = -cosc*sina + sinc*sinb*cosa;
rotmat(2,2) = cosc*cosa + sinc*sinb*sina;
rotmat(2,3) = sinc * cosb;
rotmat(3,1) = sinc*sina + cosc*sinb*cosa;
rotmat(3,2) = -sinc*cosa + cosc*sinb*sina;
rotmat(3,3) = cosc * cosb;

rotred_matrix = (redmat*rotmat)';

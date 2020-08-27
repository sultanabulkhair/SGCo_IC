
function C = cova(it,h);

%------------------------------------------
% Compute covariance for reduced distance h
%------------------------------------------
%

warning('off','all');

if (it < 1) % Nugget effect

  C = (h<eps)+0;

elseif (it < 2) % Spherical model

  C = 1 - 1.5*min(h,1) + 0.5*(min(h,1).^3);

elseif (it < 3) % Exponential model

  C = exp(-3*h);

elseif (it < 4) % Cubic model

  C = 1 - 7*(min(h,1).^2) + 35/4*(min(h,1).^3) - 7/2*(min(h,1).^5) + 3/4*(min(h,1).^7);

elseif (it < 5) % Gaussian model

  C = exp(-3*h.^2);

else

  error('Unavailable covariance model');

end

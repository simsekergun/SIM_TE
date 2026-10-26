function [n, k] = eps2nk(epsr)
% Complex relative permittivity -> refractive index n and extinction k.
%   [n, k] = eps2nk(epsr)   epsr scalar or array, exp(-i w t) convention,
%   so that epsr = (n + 1i*k)^2 with k >= 0 for passive media (imag(epsr) >= 0).
%   Example: [n,k] = eps2nk(-11.253 + 0.813i)   
% gold, 633 nm -> n = 0.1211, k = 3.3566
nt = sqrt(epsr);                 % principal branch: real(nt) >= 0
n  = real(nt);
k  = imag(nt);
if any(k(:) < 0)
    warning('eps2nk:active', 'imag(epsr) < 0 (active medium or exp(+i w t) data)');
end
end
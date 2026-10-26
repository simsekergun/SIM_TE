function X = sim_xsec(out, phi_obs, Mphi)
%SIM_XSEC  Far field, 2-D scattering width and cross sections (per unit
%          length) from the SIM boundary solution returned by sim2d.
%
%   X = sim_xsec(out)                 cross sections only
%   X = sim_xsec(out, phi_obs)        also far-field at angles phi_obs
%   X = sim_xsec(out, phi_obs, Mphi)  Mphi angles for the C_sca integral
%
%   Far-field amplitude (scattered field ~ sqrt(2/(pi kb r)) e^{i(kb r - pi/4)} T(phi)):
%     T(phi) = (i/4) \oint [ -i kb (rhat.n') U(r') - dU/dn'(r') ] e^{-i kb rhat.r'} ds'
%   Scattering width   sigma(phi) = 2 pi r |U_s|^2 = (4/kb)|T(phi)|^2
%   C_sca = (1/2pi) \int sigma dphi,   C_ext = -(4/kb) Re T(thetai)  (optical theorem)
%   C_abs = -(1/kb) \oint Im( conj(U) dU/dn ) ds     (net Poynting inflow)
%   All integrals: trapezoidal rule in t (spectrally accurate, periodic).
% 
% Ergun Simsek
% UMBC
% October, 2026

if nargin < 2
    phi_obs = [];
end
if nargin < 3
    Mphi = max(512, 2*out.N);
end
kb = out.kb;
w = (2*pi/out.N)*out.Ja;
Tfun = @(p) farfield(out, p(:).', w);
pf = (0:Mphi-1)*2*pi/Mphi;
Tf = Tfun(pf);
X.Csca = mean((4/kb)*abs(Tf).^2);
X.Cext = -(4/kb)*real(Tfun(out.thetai));
X.Cabs = -(1/kb)*sum(imag(conj(out.U).*out.dUn).*w);
X.energy_residual = (X.Cext - X.Csca - X.Cabs)/X.Cext;
if ~isempty(phi_obs)
    X.phi = phi_obs;
    X.T = Tfun(phi_obs);
    X.sigma = (4/kb)*abs(X.T).^2;
end
end

function T = farfield(out, p, w)
kb = out.kb;
r = out.r;
n = out.nrm;
c = cos(p(:));
s = sin(p(:));                       % (P x 1)
rn = c*n(1,:) + s*n(2,:);                            % rhat . n'
rr = c*r(1,:) + s*r(2,:);                            % rhat . r'
K  = (-1i*kb*rn.*repmat(out.U,numel(p),1) - repmat(out.dUn,numel(p),1)) ...
    .*exp(-1i*kb*rr);
T  = (1i/4)*(K*w.').';
end

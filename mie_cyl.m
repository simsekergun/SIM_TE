function M = mie_cyl(a, lam, epsr, phi, opts)
%MIE_CYL  Exact (Mie / Lorenz-Mie) series for a circular cylinder of radius a,
%         normally incident unit plane wave, TE_z or TM_z, complex permittivity.
%
%   M = mie_cyl(a, lam, epsr, phi)            boundary fields at angles phi
%   M = mie_cyl(a, lam, epsr, phi, opts)      opts.pol ('TE'), .thetai (0),
%                                             .epsb (1), .nmax (auto)
%
%   exp(-i w t) convention.  Incident field  U_inc = exp(i kb r cos(phi-thetai))
%                                                  = sum_n i^n J_n(kb r) e^{in(phi-thetai)}
%   Scattered      U_s  = sum_n i^n a_n H_n^(1)(kb r) e^{in(phi-thetai)}
%   Interior       U_1  = sum_n i^n c_n J_n(k1 r)     e^{in(phi-thetai)}
%
%   With m = sqrt(epsr/epsb), x = kb a, mx = k1 a:
%     TE_z (U=Hz; Hz and (1/eps)dHz/dr continuous):
%        a_n = [J_n(x) J_n'(mx) - m J_n(mx) J_n'(x)] / [m J_n(mx) H_n'(x) - H_n(x) J_n'(mx)]
%     TM_z (U=Ez; Ez and dEz/dr continuous):
%        a_n = [m J_n(x) J_n'(mx) - J_n(mx) J_n'(x)] / [J_n(mx) H_n'(x) - m H_n(x) J_n'(mx)]
%
%   Output fields: U, Us, dUn (outward d/dr, exterior side) at r = a, the
%   coefficients a_n, n, far-field amplitude handle T(phi) and the
%   cross sections per unit length Cext, Csca, Cabs (same length unit as a)
%   and efficiencies Q = C/(2a).
% 
% Ergun Simsek
% UMBC
% October, 2026

if nargin < 5
    opts = struct();
end

pol    = upper(getopt(opts,'pol','TE'));
thetai = getopt(opts,'thetai',0);
epsb   = getopt(opts,'epsb',1);

% wavenumbers;
k0 = 2*pi/lam;
kb = k0*sqrt(epsb);
% contrast
m = sqrt(epsr/epsb);
x  = kb*a;
mx = m*x;
nmax = getopt(opts,'nmax', ceil(x + 4*x^(1/3) + 50));
n  = -nmax:nmax;

% Bessel/Hankel functions and their derivatives
Jx  = besselj(n, x);
dJx  = dbess(@besselj, n, x);
Hx  = besselh(n, 1, x);
dHx  = (besselh(n-1,1,x) - besselh(n+1,1,x))/2;
Jm  = besselj(n, mx);
dJm  = dbess(@besselj, n, mx);

switch pol
    case 'TE'
        an = (Jx.*dJm - m*Jm.*dJx)./(m*Jm.*dHx - Hx.*dJm);
    case 'TM'
        an = (m*Jx.*dJm - Jm.*dJx)./(Jm.*dHx - m*Hx.*dJm);
end
cn = (Jx + an.*Hx)./Jm;                    % interior coefficients

phi = phi(:).';
E   = exp(1i*(phi.' - thetai)*n).*repmat(1i.^n, numel(phi), 1);
M.U   = (E*(Jx + an.*Hx).').';
M.Us  = (E*(an.*Hx).').';
M.dUn = kb*(E*(dJx + an.*dHx).').';
M.an = an;
M.cn = cn;
M.n = n;
M.kb = kb;
M.k1 = k0*sqrt(epsr);
M.T    = @(p) exp(1i*(p(:) - thetai)*n)*an.';          % far-field amplitude
M.Csca = (4/kb)*sum(abs(an).^2);
M.Cext = -(4/kb)*sum(real(an));
M.Cabs = M.Cext - M.Csca;
M.Qsca = M.Csca/(2*a);
M.Qext = M.Cext/(2*a);
M.Qabs = M.Cabs/(2*a);
end

function d = dbess(f, n, z)
d = (f(n-1, z) - f(n+1, z))/2;
end

function v = getopt(s, f, d)
if isfield(s,f) && ~isempty(s.(f))
    v = s.(f);
else
    v = d;
end
end

function out = sim2d(geom, lam, epsr, N, opts)
%SIM2D  Spectral integral method (SIM) for 2-D scattering of a plane wave by
%       a homogeneous, penetrable cylinder with a smooth closed boundary.
%       TE_z (H = Hz z^) or TM_z (E = Ez z^) polarization, complex permittivity.
%
%   out = sim2d(geom, lam, epsr, N)
%   out = sim2d(geom, lam, epsr, N, opts)
%
%   Time convention exp(-i*omega*t)  ->  lossy media have imag(epsr) > 0
%   (gold at 633 nm: epsr ~ -11.3 + 0.8i).  Same convention as die.m / SIM01.
%
%   INPUT
%     geom  struct describing the boundary (parameter t in [0,2*pi), CCW):
%             geom.type = 'circle'  : geom.a                (radius)
%                         'ellipse' : geom.a, geom.b        (semi-axes x,y)
%                         'pert'    : geom.a, geom.b        (shape0 type 3 of die.m)
%                         'rose'    : geom.a, geom.d, geom.m (rho=a(1+d cos(m t)))
%                         'custom'  : geom.fun = @(t) [r,dr,ddr]  (2xN each)
%     lam   vacuum wavelength (same length unit as the geometry)
%     epsr  complex relative permittivity of the cylinder (scalar)
%     N     number of boundary points (even integer)
%     opts  (optional) struct:
%             .pol    'TE' (default) or 'TM'
%             .thetai incidence angle in rad (default 0, propagation along +x)
%             .epsb   real, positive background permittivity (default 1)
%             .Lon    window switch-on threshold (default 18).  When
%                     Im(k1)*Dmax > Lon (Dmax = largest chord of the
%                     boundary) the global subtraction of SIM01/die.m would
%                     lose ~Im(k1)*Dmax/2.3 digits to cancellation, so the
%                     interior log-singularity subtraction is windowed.
%                     Inf -> never window (original SIM01 splitting).
%             .Lcut   window width parameter (default 8): subtraction is
%                     restricted to |t-t'| < Lcut/(Im(k1)*max|r'(t)|) with a
%                     C-infinity cutoff that is flat at t = t'.
%             .keepMatrices  (default false) also return out.S (Schur
%                     complement C^b - D^b (rho D^1)^{-1} C^1) and out.Dr
%                     (rho D^1) for condition-number studies.
%
%   OUTPUT  struct with fields (row vectors of length N unless noted)
%     theta  parameter values t_m          r, dr   2xN position / tangent
%     Ja     |r'(t)|                       nrm     2xN outward unit normal
%     U      total field on the boundary   (Hz for TE, Ez for TM)
%     dUn    outward normal derivative of U, exterior side
%     dUn_in outward normal derivative of U, interior side
%     Uinc   incident field on boundary    Us = U - Uinc (scattered)
%     phi,psi  Fourier coefficients (FFT ordering) of U and dUn
%     kb,k1  background / interior wavenumbers, plus input echo
%
%   The unit incident wave is U_inc = exp(i*kb*(x cos(thetai) + y sin(thetai))).
%
%   Boundary conditions used (see accompanying LaTeX note):
%     TM_z :  Ez continuous,  dEz/dn continuous             (mu_r = 1)
%     TE_z :  Hz continuous,  (1/eps) dHz/dn continuous
%             -> dHz/dn|in = (epsr/epsb) dHz/dn|out
%
%   Based on J. Liu and Q. H. Liu, "A spectral integral method for periodic
%   and nonperiodic structures" 2004.
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
Lcut   = getopt(opts,'Lcut',8);
Lon    = getopt(opts,'Lon',18);
if mod(N,2)~=0 
    error('sim2d:N','N must be even'); 
end
if ~(isreal(epsb) && epsb>0) 
    error('sim2d:epsb','background must be lossless'); 
end

% wavenumbers
k0 = 2*pi/lam;
kb = k0*sqrt(epsb);
k1 = k0*sqrt(epsr);            % principal branch: imag(k1) >= 0 for imag(epsr) >= 0
switch pol
    case 'TE', rho = epsr/epsb;  % dU/dn|in = rho * dU/dn|out
    case 'TM', rho = 1;
    otherwise, error('sim2d:pol','pol must be TE or TM');
end

%---- discretisation ----
ii  = 1:N;
th  = (ii-1)*2*pi/N;                 % basis = testing points (t_m)
[rr,dr,ddr] = curve(geom, th);
Ja  = sqrt(dr(1,:).^2 + dr(2,:).^2);
kap = (dr(1,:).*ddr(2,:) - dr(2,:).*ddr(1,:))./Ja.^3;   % signed curvature
nrm = [dr(2,:); -dr(1,:)]./[Ja; Ja];                      % outward normal (CCW)

% window half-widths in t for exterior (index 1) and interior (index 2)
tw = [Inf, Inf];
Rmax = max(max(hypot(rr(1,:).' - rr(1,:), rr(2,:).' - rr(2,:))));
if isfinite(Lon) && imag(k1)*Rmax > Lon
    tw(2) = Lcut/(imag(k1)*max(Ja));
end

%---- incident field ----
inc = exp(1i*kb*(cos(thetai)*rr(1,:) + sin(thetai)*rr(2,:)));
ff  = (fft(inc, N)/N).';              % Fourier coefficients, inc(t) = sum_k ff_k e^{+ikt}
                                      % (same e^{+ikt} convention as the matrix rows; the
                                      %  former ifft(inc) used e^{-ikt}, which is exact only for
                                      %  boundaries symmetric under t -> -t)

%---- Fourier coefficients of ln|2 sin(t/2)| : a_0 = 0, a_n = -1/(2|n|) --
mf = 2*N;  
nn = 1:N-1;
a0 = -1.0./(2.0*ii);

A = zeros(N); 
B = A; 
C = A; 
D = A; 
EE = A;
kk = [kb, k1];
for j = 1:N                               % loop over testing points
    tm  = th(j);
    etm = exp(1i*(ii-1)*tm);              % e^{i n t_m}, FFT ordering
    rt  = rr(:,j);

    aj = zeros(1,mf);                     % zero-padded log coefficients
    aj(2:N)    = a0(nn).*exp(1i*nn*tm);
    aj(N+1:mf) = a0(N:-1:1).*exp(-1i*(N:-1:1)*tm);
    ajk = fft(aj, mf);

    for s = 1:2                           % s=1 exterior, s=2 interior
        k = kk(s);
        [g,h,p0,q0] = kernel_row(k, th, rr, Ja, nrm, kap, tm, rt, tw(s), j);
        u = padconv(p0, ajk, N);
        v = padconv(q0, ajk, N);
        hk = -2*(2*pi*h + k.*v);          % matrix row of the double layer
        gk =  2*(-u + 2*pi*g);            % matrix row of the single layer
        if s == 1
            A(j,:) = hk + etm;  
            B(j,:) = gk;
        else
            C(j,:) = hk - etm;  
            D(j,:) = gk;
        end
    end
    EE(j,:) = etm;
end

%---- solve  [A B; C rho*D] [phi; psi] = [2F; 0] ------------------------
Dr  = rho*D;
X   = Dr\C;
phi = (A - B*X)\(EE*(2*ff));
psi = -X*phi;

U   = N*ifft(phi.', N);               % U(t_m) = sum_k phi_k e^{+ik t_m}
dUn = N*ifft(psi.', N);

out = struct('theta',th,'r',rr,'dr',dr,'Ja',Ja,'kappa',kap,'nrm',nrm, ...
    'U',U,'dUn',dUn,'dUn_in',rho*dUn,'Uinc',inc,'Us',U-inc, ...
    'phi',phi.','psi',psi.','pol',pol,'lam',lam,'epsr',epsr,'epsb',epsb, ...
    'k0',k0,'kb',kb,'k1',k1,'thetai',thetai,'N',N,'tw',tw,'geom',geom);
if getopt(opts,'keepMatrices',false)       % optional: system matrices for
    out.S  = A - B*X;                      % condition-number studies
    out.Dr = Dr;                           % (Schur complement, rho*D^1)
end
end

%=========================================================================
function [g,h,p,q] = kernel_row(k, th, rr, Ja, nrm, kap, tm, rt, tw, jself)
% Fourier coefficients (over source parameter t') of the smooth parts of
% G and dG/dn' times |r'(t')| and of the two functions multiplying ln|2sin|.
%   G      = (i/4) H0(kR)
%   dG/dn' = -(i k/4) H1(kR) (r'-r).n'/R
%   Gbar   = G      + (1/2pi)   ln|2 sin((t-t')/2)| J0(kR) chi
%   Hbar   = dG/dn' - (k/2pi)   ln|2 sin((t-t')/2)| J1(kR) (r'-r).n'/R chi
% chi = 1 (no window) or a C-infinity cutoff that is identically 1 to all
% orders at t=t' (so Gbar, Hbar remain C-infinity).
euler = 0.57721566490153288;
d    = th - tm;
self = false(size(th)); 
self(jself) = true;
Rx = rr(1,:) - rt(1);  
Ry = rr(2,:) - rt(2);
R  = sqrt(Rx.^2 + Ry.^2);  
R(self) = 1;            % dummy, overwritten
RnR = (Rx.*nrm(1,:) + Ry.*nrm(2,:))./R;            % (r'-r).n'/R (signed)
x  = k*R;

xl = zeros(size(th));
xl(~self) = log(2*abs(sin(d(~self)/2)));

chi = cutoff(d, tw);
act = chi > 0;
J0 = zeros(size(th));  
J1 = J0;
J0(act) = besselj(0, x(act));
J1(act) = besselj(1, x(act));

H0 = besselh(0, 1, x);
H1 = besselh(1, 1, x);

p = J0.*chi;          
p(self) = 1;
q = J1.*RnR.*chi;     
q(self) = 0;

g = (1i/4)*H0 + (1/(2*pi))*xl.*p;
g(self) = -(euler + log(k*Ja(self)/2))/(2*pi) + 1i/4;

h = -(1i/4)*k*RnR.*H1 - (k/(2*pi))*xl.*q;
h(self) = -kap(self)/(4*pi);                        % limit R -> 0

N = numel(th);
g = ifft(g.*Ja, N);  
h = ifft(h.*Ja, N);
p = ifft(p.*Ja, N);  
q = ifft(q.*Ja, N);
end

%=========================================================================
function u0 = padconv(p0, ajk, N)
% discrete convolution (in Fourier index) of p0 with the log coefficients
mf = 2*N;
p = zeros(1,mf);
p(1:N/2) = p0(1:N/2);
p(mf-N/2+1:mf)  = p0(N/2+1:N);
u = ifft(fft(p,mf).*ajk, mf);
u0 = [u(1:N/2), u(mf-N/2+1:mf)];
end

%=========================================================================
function c = cutoff(t, tw)
% C-infinity cutoff: 1 at t=0 (flat to all orders), 0 for |t| >= tw.
if ~isfinite(tw), c = ones(size(t)); return; end
t  = abs(angle(exp(1i*t)));                % wrap to [0, pi]
s  = min(t/tw, 1);
f1 = zeros(size(s));  
f2 = f1;
m1 = s > 0;  f1(m1) = exp(-1./s(m1));
m2 = s < 1;  f2(m2) = exp(-1./(1-s(m2)));
c  = f2./(f1 + f2);
end

%=========================================================================
function [r,dr,ddr] = curve(geom, t)
switch lower(geom.type)
    case 'circle'
        a = geom.a;
        r   = [ a*cos(t);  a*sin(t)];
        dr  = [-a*sin(t);  a*cos(t)];
        ddr = [-a*cos(t); -a*sin(t)];
    case 'ellipse'
        a = geom.a; 
        b = geom.b;
        r   = [ a*cos(t);  b*sin(t)];
        dr  = [-a*sin(t);  b*cos(t)];
        ddr = [-a*cos(t); -b*sin(t)];
    case 'pert'                               % shape0(3,...) of die.m
        a = geom.a; 
        b = geom.b;
        s = t + 0.05*sin(2*t);  
        f = 1 + 0.1*cos(2*t);
        r   = [ a*cos(s);  b*sin(s)];
        dr  = [-a*sin(s).*f;  b*cos(s).*f];
        ddr = [-r(1,:).*f.^2 + 0.2*a*sin(s).*sin(2*t); ...
               -r(2,:).*f.^2 - 0.2*b*cos(s).*sin(2*t)];
    case 'rose'                               % rho = a(1 + d cos(m t))
        a = geom.a; 
        dd = geom.d; 
        m = geom.m;
        rho  = a*(1 + dd*cos(m*t));
        drho = -a*dd*m*sin(m*t);
        d2   = -a*dd*m^2*cos(m*t);
        r   = [rho.*cos(t); rho.*sin(t)];
        dr  = [drho.*cos(t) - rho.*sin(t); drho.*sin(t) + rho.*cos(t)];
        ddr = [(d2 - rho).*cos(t) - 2*drho.*sin(t); ...
               (d2 - rho).*sin(t) + 2*drho.*cos(t)];
    case 'custom'
        [r,dr,ddr] = geom.fun(t);
    otherwise
        error('sim2d:geom','unknown geometry type %s', geom.type);
end
end

function v = getopt(s, f, d)
    if isfield(s,f) && ~isempty(s.(f)) 
        v = s.(f); 
    else 
        v = d; 
    end
end

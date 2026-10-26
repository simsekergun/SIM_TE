function e = eps_gold(lam_nm)
%EPS_GOLD  Relative permittivity of gold, 
% exp(-i w t) convention (Im e > 0).
%
%   e = eps_gold(lam_nm)
%       Analytic Drude + two critical-point model fitted to Johnson & Christy:
%       P. G. Etchegoin, E. C. Le Ru, M. Meyer, J. Chem. Phys. 125, 164705
%       (2006); erratum 127, 189901 (2007).  Valid roughly 400-1000 nm.

lam_nm = double(lam_nm);
einf = 1.53;  
lp = 145;  
gp = 17000;                  % nm
A  = [0.94, 1.36];   
ph = [-pi/4, -pi/4];
li = [468, 331];     
gi = [2300, 940];                % nm
l  = lam_nm;
e  = einf - 1./(lp^2*(1./l.^2 + 1i./(gp*l)));
for p = 1:2
    e = e + A(p)/li(p)*( exp(1i*ph(p))./(1/li(p) - 1./l - 1i/gi(p)) ...
                       + exp(-1i*ph(p))./(1/li(p) + 1./l + 1i/gi(p)) );
end
end

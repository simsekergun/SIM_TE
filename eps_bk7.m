function e = eps_bk7(lam_nm)
%EPS_BK7  Relative permittivity of N-BK7 glass 
% (Schott Sellmeier formula, lossless in the visible), 
% lam in nm.  e = n^2.
L2 = (double(lam_nm)/1000).^2;                        % um^2
B = [1.03961212, 0.231792344, 1.01046945];
C = [0.00600069867, 0.0200179144, 103.560653];        % um^2
e = 1 + B(1)*L2./(L2-C(1)) + B(2)*L2./(L2-C(2)) + B(3)*L2./(L2-C(3));
end

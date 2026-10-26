%VALIDATE_TE  Validation of the TE_z (and TM_z) spectral integral method
%             against the exact Mie series for glass (BK7) and gold cylinders.
%   Run from this folder:  >> validate_te
%   Lengths in nm.  Expected results (from the NumPy port of this code) are
%   quoted in the comments; MATLAB should reproduce them to ~1e-13.
% 
% Ergun Simsek
% UMBC
% October, 2026

clear; 

%% 1. Surface-field convergence vs N (TE and TM)
cases = { 'Glass BK7, a=1000 nm, 633 nm', 1000, 633, eps_bk7(633);
          'Gold, a=50 nm, 520 nm',          50,  520, eps_gold(520);
          'Gold, a=100 nm, 633 nm',        100,  633, eps_gold(633);
          'Gold, a=200 nm, 800 nm',        200,  800, eps_gold(800)};
Ns = [16 24 32 48 64 96 128 192 256 384 512];
err = zeros(size(cases,1), numel(Ns), 2);
pols = {'TE','TM'};
for c = 1:size(cases,1)
    a = cases{c,2}; 
    lam = cases{c,3}; 
    ep = cases{c,4};
    for ip = 1:2
        for k = 1:numel(Ns)
            o = sim2d(struct('type','circle','a',a), lam, ep, Ns(k), struct('pol',pols{ip}));
            M = mie_cyl(a, lam, ep, o.theta, struct('pol',pols{ip}));
            err(c,k,ip) = max(abs(o.U - M.U))/max(abs(M.U));
        end
        fprintf('%-30s %s: ', cases{c,1}, pols{ip}); fprintf('%8.1e', err(c,:,ip)); fprintf('\n');
    end
end
figure(1); clf;
semilogy(Ns, err(:,:,1).', 'o-'); 
grid on;
xlabel('{\it{N}}'); 
ylabel('Relative Error'); 
legend(cases(:,1), 'Location','east');

% window off (original SIM01 splitting) for the large gold cylinder
e0 = zeros(size(Ns));
for k = 1:numel(Ns)
    o = sim2d(struct('type','circle','a',500), 800, eps_gold(800), Ns(k), struct('Lon',Inf));
    M = mie_cyl(500, 800, eps_gold(800), o.theta);
    e0(k) = max(abs(o.U - M.U))/max(abs(M.U));
end
fprintf('Gold a=500 nm, no window      TE: '); 
fprintf('%8.1e', e0); 
fprintf('\n');   % O(1) errors
hold on; 
semilogy(Ns, e0, 'x--');

%% 2. Surface field and far field at 633 nm (TE)
geo = {struct('type','circle','a',200), struct('type','circle','a',200)};
eps633 = [eps_bk7(633), eps_gold(633)];  
aa = [200 200];  
NN = [64 96];
phis = linspace(0, 2*pi, 721);
figure(2); clf;
for c = 1:2
    o = sim2d(geo{c}, 633, eps633(c), NN(c));
    th = linspace(0,2*pi,721);  
    M = mie_cyl(aa(c), 633, eps633(c), th);
    subplot(2,2,c); 
    plot(th*180/pi, abs(M.U), 'k-', o.theta*180/pi, abs(o.U), 'ro');
    xlabel('\theta (deg)'); 
    ylabel('|{\it{H}}_z|'); 
    legend('Mie','SIM','Location','South'); 
    xlim([0 360]);
    X  = sim_xsec(o, phis);
    sm = (4/M.kb)*abs(M.T(phis)).^2;
    grid on;
    subplot(2,2,c+2); 
    semilogy(phis*180/pi, sm, 'k-', phis(1:12:end)*180/pi, X.sigma(1:12:end), 'ro');
    xlabel('\theta (deg)'); 
    ylabel('\sigma (nm)'); 
    grid on;
    xlim([0 360]);
    fprintf('far field case %d: max rel err %.1e\n', c, max(abs(X.sigma - sm.'))/max(sm));
end

%% 3. Cross-section spectra 400-900 nm (N = 64)
lams = 400:10:800;
mats = {'gold', 200, @eps_gold; 'glass', 200, @eps_bk7};
fig = figure(3); clf;
fig.Position = [100 100 1000 400];
for c = 1:2
    a = mats{c,2}; 
    f = mats{c,3};
    Qs = zeros(numel(lams),3); 
    Qm = Qs; 
    res = zeros(numel(lams),1);
    for k = 1:numel(lams)
        o = sim2d(struct('type','circle','a',a), lams(k), f(lams(k)), 64);
        X = sim_xsec(o);  
        M = mie_cyl(a, lams(k), f(lams(k)), 0);
        Qs(k,:) = [X.Cext X.Csca X.Cabs]/(2*a);  
        Qm(k,:) = [M.Qext M.Qsca M.Qabs];
        res(k) = X.energy_residual;
    end
    fprintf('%s TE spectra: max|dQext|=%.1e max|dQabs|=%.1e max|energy res|=%.1e\n', ...
        mats{c,1}, max(abs(Qs(:,1)-Qm(:,1))), max(abs(Qs(:,3)-Qm(:,3))), max(abs(res)));
    subplot(1,2,c); 
    plot(lams, Qm(:,1), 'k-', lams, Qs(:,1), 'ro', lams, Qm(:,3), 'b-', lams, Qs(:,3), 'bs');
    xlabel('\lambda (nm)'); 
    ylabel('{\it{Q}}'); 
    grid on;
    legend('{\it{Q}}_{ext} Mie','{\it{Q}}_{ext} SIM','{\it{Q}}_{abs} Mie','{\it{Q}}_{abs} SIM');
end

%% 4. Non-convex gold "rose" (no analytic solution): self-convergence + optical theorem
rose = struct('type','rose','a',150,'d',0.25,'m',5);
ref = sim2d(rose, 633, eps_gold(633), 768);
for N = [64 96 128 192 256 384]
    o = sim2d(rose, 633, eps_gold(633), N);  X = sim_xsec(o);
    e = max(abs(o.U - ref.U(1:768/N:end)))/max(abs(ref.U));
    fprintf('rose N=%3d  self-err %.1e  (Cext-Csca-Cabs)/Cext = %.1e\n', N, e, X.energy_residual);
end
% expected: self-err 4e-13 and energy residual 4e-14 at N = 384

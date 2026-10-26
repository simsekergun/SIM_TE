clear;
set(0, 'DefaultFigureCreateFcn', @(h,~) set(h, 'Renderer', 'painters'))

geo = {struct('type','circle','a',1000), struct('type','circle','a',100)};
eps633 = [eps_gold(633)];  
aa = [1000 100];  
NN = [256 512];
phis = linspace(0, 2*pi, 721);
figure(2); clf;
for c = 1
    o = sim2d(geo{c}, 633, eps633(c), NN(c));
    th = linspace(0,2*pi,721);  
    M = mie_cyl(aa(c), 633, eps633(c), th);
    subplot(2,1,c); plot(th*180/pi, abs(M.U), 'k-', o.theta*180/pi, abs(o.U), 'ro');
    xlabel('\theta (deg)'); 
    ylabel('|H_z|'); 
    legend('Mie','SIM','Location','South'); 
    xlim([0 360]);
    X  = sim_xsec(o, phis);
    sm = (4/M.kb)*abs(M.T(phis)).^2;
    grid on;
    subplot(2,1,c+1); 
    semilogy(phis*180/pi, sm, 'k-', phis(1:12:end)*180/pi, X.sigma(1:12:end), 'ro');
    xlabel('\phi (deg)'); 
    ylabel('\sigma (nm)'); 
    xlim([0 360]);
    grid on;
    print -dpng figure_Hz_vs_angle
    fprintf('far field case %d: max rel err %.1e\n', c, max(abs(X.sigma - sm.'))/max(sm));
end

%% 3. Cross-section spectra 400-900 nm (N = 64)
lams = 400:10:900;
mats = {'gold', 50, @eps_gold; 'glass', 500, @eps_bk7};
figure(3); clf;
for c = 1
    a = mats{c,2}; f = mats{c,3};
    Qs = zeros(numel(lams),3); 
    Qm = Qs; 
    res = zeros(numel(lams),1);
    for k = 1:numel(lams)
        o = sim2d(struct('type','circle','a',a), lams(k), f(lams(k)), 64);
        X = sim_xsec(o);  M = mie_cyl(a, lams(k), f(lams(k)), 0);
        Qs(k,:) = [X.Cext X.Csca X.Cabs]/(2*a);  
        Qm(k,:) = [M.Qext M.Qsca M.Qabs];
        res(k) = X.energy_residual;
    end
    fprintf('%s TE spectra: max|dQext|=%.1e max|dQabs|=%.1e max|energy res|=%.1e\n', ...
        mats{c,1}, max(abs(Qs(:,1)-Qm(:,1))), max(abs(Qs(:,3)-Qm(:,3))), max(abs(res)));
    plot(lams, Qm(:,1), 'k-', lams, Qs(:,1), 'ro', lams, Qm(:,3), 'b-', lams, Qs(:,3), 'bs');
    xlabel('\lambda (nm)'); 
    ylabel('Q'); 
    legend('Q_{ext} Mie','Q_{ext} SIM','Q_{abs} Mie','Q_{abs} SIM');
    grid on;
    print -dpng figure_Q_vs_lambda
end
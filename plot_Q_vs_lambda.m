clear

%% 3. Cross-section spectra 400-900 nm (N = 64)
lams = 400:10:800;
mats = {'gold', 200, @eps_gold};
fig = figure(3); clf;
fig.Position = [200 300 700 500];
c = 1;
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
    plot(lams, Qm(:,1), 'k-', lams, Qs(:,1), 'ro', lams, Qm(:,3), 'b-', lams, Qs(:,3), 'bs');
    xlabel('\lambda (nm)'); 
    ylabel('{\it{Q}}'); 
    grid on;
    legend('{\it{Q}}_{ext} Mie','{\it{Q}}_{ext} SIM','{\it{Q}}_{abs} Mie','{\it{Q}}_{abs} SIM');
  %  title(sprintf('%s, a = %d nm, TE_z', mats{c,1}, a));
  print -dpng figure_Q_vs_lambda

clear; % close all;
set(0, 'DefaultFigureCreateFcn', @(h,~) set(h, 'Renderer', 'painters'))
%% 1. Surface-field convergence vs N (TE and TM)
cases = { 'Au, Rose, r = 200 nm, \lambda = 633 nm', 200, 100, 2, 633, eps_gold(633)};
% Ns = [ 16 20 24 28 32 48 64 72 96 128];
Ns = [ 16 20 24 ];
err = zeros(size(cases,1), numel(Ns), 2);
pols = {'TE','TM'};

for c = 1:size(cases,1)
    a = cases{c,2}; 
    lam = cases{c,3}; 
    ep = cases{c,4};
    n = sqrt((sqrt(real(ep).^2+imag(ep).^2)+real(ep))/2);
    circumf = 2*pi*a/lam;
    for ip = 1:2
        for k = 1:numel(Ns)
            o = sim2d(struct('type','circle','a',a), lam, ep, Ns(k), struct('pol',pols{ip}));
            M = mie_cyl(a, lam, ep, o.theta, struct('pol',pols{ip}));
            err(c,k,ip) = max(abs(o.U - M.U))/max(abs(M.U));
        end
        fprintf('%-30s %s: ', cases{c,1}, pols{ip}); fprintf('%8.1e', err(c,:,ip)); fprintf('\n');
    end
end
figure(2); clf;
semilogy(Ns/circumf, err(:,:,1).', 'o-'); grid on;
xlabel('Points per wavelength'); 
ylabel('Maximum Relative Error'); 
grid on;
%title('TE_z convergence');
legend(cases(:,1), 'Location','northeast');
print -dpng figure_error_vs_ppw_au

ppw = Ns/circumf;
errors = err(:,:,1).';
save Au_results 'ppw' 'errors'
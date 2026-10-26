clear;

load die_results 
die_ppw = ppw;
die_errors =  errors;

load au_results 
au_ppw = ppw;
au_errors =  errors;

figure(9); clf;
semilogy(die_ppw, die_errors,'-s', au_ppw, au_errors,'-d');
xlabel('Points per wavelength');
ylabel('Relative Error');
grid on;
legend('BK7 glass', 'Au')
print -dpng figure_error_vs_ppw
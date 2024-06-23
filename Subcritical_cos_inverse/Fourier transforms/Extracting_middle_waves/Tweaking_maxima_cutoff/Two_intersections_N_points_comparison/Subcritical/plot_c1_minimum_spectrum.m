[~,min_i]=min(nonlinear_c1_power);

%min_i=min_i-1

figure; clf; hold on;
stem(nonlinear_freqs{min_i},nonlinear_yshift{min_i})
xline(ks(min_i),'--k')
xlim([0,5*ks(min_i)])
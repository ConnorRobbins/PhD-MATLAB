norm1=load('lowamp_support15_for_comparison.mat').Theta_b_norm;
norm2=load('wavenumber1.mat').Theta_b_norm;



figure(2); clf; hold on;
plot(log10(norm1),'-r')
plot(log10(norm2),'-b')
legend('comparison','single peak')
ylabel('log10(norm(Theta b))')
xlabel('rank')
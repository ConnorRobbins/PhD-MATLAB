
load('unforced_inverse_max_data.mat','Ns','maxvalstore')


savefigure=1;
%%



lineN=linspace(2.2,3.3,100);
yoffset=2.1;
lineY=-2*lineN +yoffset;


figure(4680); clf; hold on; box on;
plot(lineN,lineY,':k',LineWidth=1.5)
plot(log10(Ns),log10(maxvalstore),'-xr',LineWidth=1.5,MarkerSize=10)
ylabel('$\log_{10}(\max(|y_b|))$',Interpreter='latex',FontSize=18)
xlabel('$\log_{10}(N)$',Interpreter='latex',FontSize=18)
%legend(Interpreter="latex",FontSize=14,Location='southwest')


if savefigure==1
    saveas(4680,"unforced_solitary_convergence.eps",'epsc')
    saveas(4680,"unforced_solitary_convergence",'fig')
end
s1=load("a1E-3.mat");
s2=load("a1E-2.mat");
s3=load("a2E-2.mat");





Fplot=linspace(1,1.2,1000);
unforced=1+2*(Fplot-1);


figure(1); clf; hold on; box on;

%plot(Fplot,unforced,'--k',LineWidth=2,HandleVisibility='off')
%plot(Fplot,unforced,'--k',LineWidth=2,DisplayName="Unforced solitary"+ newline + "wave")

plot(s1.F_lower,s1.max_y_lower,'-b',LineWidth=2,DisplayName="$a=1\times 10^{-3}$")
plot(s1.F_upper,s1.max_y_upper,'-b',LineWidth=2,HandleVisibility='off')

plot(s2.F_lower,s2.max_y_lower,'-r',LineWidth=2,DisplayName="$a=1\times 10^{-2}$")
plot(s2.F_upper,s2.max_y_upper,'-r',LineWidth=2,HandleVisibility='off')

plot(s3.F_lower,s3.max_y_lower,'-k',LineWidth=2,DisplayName="$a=2\times 10^{-2}$")
plot(s3.F_upper,s3.max_y_upper,'-k',LineWidth=2,HandleVisibility='off')




xlabel('$F$',Interpreter='latex',FontSize=18)
ylabel('$y_f(0)$',Interpreter='latex',FontSize=18)
legend(Location='east',Interpreter="latex",FontSize=14)
box on
%xlim([Phi_t(1),Phi_t(end)])
ylim([1,1.25])
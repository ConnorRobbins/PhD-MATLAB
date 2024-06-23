N=41;
L=20;
F=0.8;
a=-0.01;
b=1;

ijac=0; %supercritical doesn't seem to like recalculating jacobian

x=linspace(-L,L,N);
del_x=2*L/(N-1);


forcing=a*exp(-((b*x).^2));
%forcing=a*exp(-((b*(x-5)).^2)) +a*exp(-((b*(x+5)).^2));


eta_init_guess=zeros(1,N);



mu=F-1;






eta=FUNCTION_kdv_newton_finite_difference_jac(N,forcing,mu,del_x,eta_init_guess,ijac);

figure(1); clf; hold on;
plot(x,eta,'-b')
%plot(x,eta_init,'-r')


figure(2); clf; hold on;
plot(x,1+eta,'-b')
plot(x,forcing,'-r')




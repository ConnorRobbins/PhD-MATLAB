clc

N=541;
L=20;
F=1.2;
a=0.01;
b=1;

ijac=1; %supercritical doesn't seem to like recalculating jacobian

x=linspace(-L,L,N);
del_x=2*L/(N-1);


forcing=a*exp(-((b*x).^2));
%forcing=a*exp(-((b*(x-5)).^2)) +a*exp(-((b*(x+5)).^2));
%forcing= 0.5*a*(tanh(10-x)+tanh(10+x));
% forcing= 0.5*a*(tanh(10-x)+tanh(10+x))-0.5*a*exp(-((b*x).^2));
% forcing=0.5*a*(tanh(10-x)+tanh(10+x)).*cos(5*x);


eta_init_guess=zeros(1,N);




mu=F-1;






eta=FUNCTION_kdv_newton_explicit_jac(N,forcing,mu,del_x,eta_init_guess,ijac);

figure(1); clf; hold on;
plot(x,eta,'-r')
%plot(x,eta_init,'-r')


figure(2); clf; hold on;
plot(x,1+eta,'-b')
plot(x,forcing,'-r')

figure(3); clf; hold on;
subplot(2,1,1)
plot(x,eta,'-b')
ylabel('eta')
subplot(2,1,2)
plot(x,forcing,'-r')
ylabel('forcing')




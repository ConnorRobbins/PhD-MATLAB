clc

N=241;
L=30;
F=0.8;
ijac=1; 


%%%% forcing of form f(x)= a*exp(-((b*(x-c)).^2)) + a*exp(-((b*(x+c)).^2))
a=-0.04;
b=1;
c_init=5;

%%%% initial guess for surface
eta_init=zeros(1,N);




x=linspace(-L,L,N);
del_x=2*L/(N-1);
mu=F-1;
initial_unknowns=[eta_init,c_init];





[eta,c]=FUNCTION_kdv_subcritical_trapped_waves(N,a,b,mu,x,del_x,ijac,initial_unknowns);
forcing=a*exp(-((b*(x-c)).^2)) +a*exp(-((b*(x+c)).^2));

figure(1); clf; hold on;
plot(x,eta,'-r')
%plot(x,eta_init,'-r')


figure(2); clf; hold on;
plot(x,1+eta,'-b')
plot(x,forcing,'-r')




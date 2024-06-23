clc

N=321;
L=20;
F=0.8;
a=-0.01;
b=1;

ijac=1;

x=linspace(-L,L,N);
del_x=2*L/(N-1);


forcing=a*exp(-((b*x).^2));
%forcing=a*exp(-((b*(x-5)).^2)) +a*exp(-((b*(x+5)).^2));
%forcing=0.5*a*(tanh(x+5)-tanh(x-5));
%forcing=0.5*a*(tanh(x+5)-tanh(x-5)).*cos(2*x);



eta=FUNCTION_initial_guess_from_kdv(N,forcing,F,del_x,ijac);

figure(1); clf; hold on;
plot(x,eta,'-r')
%plot(x,eta_init,'-r')


figure(2); clf; hold on;
plot(x,1+eta,'-b')
plot(x,forcing,'-r')




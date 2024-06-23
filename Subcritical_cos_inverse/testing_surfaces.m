N=671;
L=25;



amplitude=0.01;
wave_k=5;
surface_support=15;
support_gradient=1;
tanh_offset=0.5;




Phi_s=linspace(-L,L,N);

tanh_support=0.5*(tanh(support_gradient*(surface_support-tanh_offset-Phi_s))+tanh(support_gradient*(surface_support-tanh_offset+Phi_s)));

Ys=1+amplitude*cos(wave_k*pi*Phi_s/(2*surface_support)).*tanh_support;



figure(2); clf;
plot(Phi_s,tanh_support)


figure(1); hold on;
plot(Phi_s,Ys)
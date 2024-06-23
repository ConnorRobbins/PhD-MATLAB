N=671;
L=25;

ks_number=4;

amplitude=0.01;
%wave_ks=linspace(1,10,ks_number);
wave_ks=1+4.*(0:1:ks_number-1);
surface_support=15;
support_gradient=2;
tanh_offset=0.5;



Phi_s=linspace(-L,L,N);
Ys=zeros(N,ks_number);
tanh_support=0.5*(tanh(support_gradient*(surface_support-tanh_offset-Phi_s))+tanh(support_gradient*(surface_support-tanh_offset+Phi_s)));


for i=1:ks_number
    wave_k=wave_ks(i);
    Ys(:,i)=1+amplitude*cos(wave_k*pi*Phi_s/(2*surface_support)).*tanh_support;
end


figure(2); clf;
plot(Phi_s,tanh_support)


figure(1); hold on;
plot(Phi_s,Ys)
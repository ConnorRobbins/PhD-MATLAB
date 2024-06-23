clear

A=0.01;
wave_k=1;

length_support_multiplier=2.5;
support_gradient=1;

%L=9.5;
N=601;
F=0.8;



surface_support=11*pi/(2*wave_k);
L=length_support_multiplier*surface_support;
Phi_s=linspace(-L,L,N);
tanh_support=0.5*(tanh(support_gradient*(surface_support-Phi_s))+tanh(support_gradient*(surface_support+Phi_s)));
eta=A*cos(wave_k*Phi_s).*tanh_support;

del_x=Phi_s(2)-Phi_s(1);
eta_xx=[0,eta(3:end)-2*eta(2:end-1)+eta(1:end-2),0]/del_x^2;

kdv_Yt=2*(F-1)*eta - (eta_xx/3)  - 1.5*eta.^2;





[frequencies,yshift] = fun_FT_topography(kdv_Yt,L,N);
[pos_frequencies,pos_yshift] = fun_FT_post_clean(frequencies,yshift);




figure(1); clf; hold on; 
plot(Phi_s,kdv_Yt)


figure(2); clf; hold on;
stem(pos_frequencies,pos_yshift)
xlim([0,5*wave_k])

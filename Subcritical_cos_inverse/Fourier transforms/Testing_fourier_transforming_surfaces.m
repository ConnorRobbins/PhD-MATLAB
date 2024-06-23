N=2070;
amplitude=0.01;


support_gradient=1;
length_support_multiplier=2.5;
wave_ks=[1,2,3];



wave_k_n=numel(wave_ks);

for i =1:wave_k_n
    wave_k=wave_ks(i);
    surface_support=9*pi/(2*wave_k);
    L=length_support_multiplier*surface_support;
    Ls(i)=L;
    Phi_s=linspace(-L,L,N);
    tanh_support=0.5*(tanh(support_gradient*(surface_support-Phi_s))+tanh(support_gradient*(surface_support+Phi_s)));
    Ys=1+amplitude*cos(wave_k*Phi_s).*tanh_support;

    [frequencies,yshift] = fun_FT_topography(Ys-1,L,N);


    figure(2*i);  clf;
    plot(Phi_s,Ys)
    title(num2str(wave_k))



    figure(2*i +1); clf;
    %stem(frequencies,yshift)
    plot(frequencies,yshift)
    title(num2str(wave_k))
end


figure(22); clf;
plot(Phi_s,tanh_support)


figure(23); clf; hold on;
plot(Phi_s,Ys)
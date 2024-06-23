clear

A=0.1;
k=1.6;


L=19.5;
N=601;
F=1.2;


x=linspace(-L,L,N);

%%%%% Assuming a forcing Acos(kx)
kdv_Yt=2*(F-1)*A*cos(k*x)   +   k^2*A*(1/3)*cos(k*x)   -   (3*A^2*((cos(k*x)).^2)/2);
%kdv_Yt=cos(k*x)+cos(2*k*x);






[frequencies,yshift] = fun_FT_topography(kdv_Yt,L,N);
[pos_frequencies,pos_yshift] = fun_FT_post_clean(frequencies,yshift);




figure(1); clf; hold on; 
plot(x,kdv_Yt)


figure(2);  clf; %hold on;
stem(pos_frequencies,pos_yshift)
xlim([0,5*k])

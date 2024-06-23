clear

A=0.1;
k=1;


L=40.;
N=1851;
F=0.8;


x=linspace(-L,L,N);


c_k0=-3*(A^2)/4;
c_k1=A*(2*(F-1) + (k^2/3));
c_k2=-3*(A^2)/4;

%%%%% Assuming a forcing Acos(kx)
kdv_Yt=c_k0+c_k1*cos(k*x)   +   c_k2*cos(2*k*x);



[frequencies,yshift] = fun_FT_topography(kdv_Yt,L,N);
[pos_frequencies,pos_yshift] = fun_FT_post_clean(frequencies,yshift);




figure(1); clf; hold on; 
plot(x,kdv_Yt)


figure(2); clf; hold on;
stem(pos_frequencies,pos_yshift)
xlim([0,5*k])



%checks
stem(0,abs(c_k0),'--r')
stem(k,abs(c_k1),'--r')
stem(2*k,abs(c_k2),'--r')
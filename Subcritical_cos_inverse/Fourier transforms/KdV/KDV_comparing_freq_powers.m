k=linspace(0.0001,3,10000);
A=0.2;
F=0.95;


c_k1=abs(A*(2*(F-1)+((k.^2)/3)));
c_k0=abs((-3*A^2)/4);

figure(10); clf; hold on; 
plot(k,c_k1)
yline(c_k0,'--k')
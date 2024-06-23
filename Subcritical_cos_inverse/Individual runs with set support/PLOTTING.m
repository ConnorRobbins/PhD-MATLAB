data1=load('support5.mat');
data2=load('support10.mat');
data3=load('support15.mat');
data4=load('support20.mat');

Phi_s=data1.Phi_s;

Ys1=data1.Ys;
Yt1=data1.Yt;
Ys2=data2.Ys;
Yt2=data2.Yt;
Ys3=data3.Ys;
Yt3=data3.Yt;
Ys4=data4.Ys;
Yt4=data4.Yt;



figure(1); clf;
subplot(2,2,1); hold on;
plot(Phi_s,Ys1,'-b')
plot(Phi_s,Yt1,'-r')
subplot(2,2,2); hold on;
plot(Phi_s,Ys2,'-b')
plot(Phi_s,Yt2,'-r')
subplot(2,2,3); hold on;
plot(Phi_s,Ys3,'-b')
plot(Phi_s,Yt3,'-r')
subplot(2,2,4); hold on;
plot(Phi_s,Ys4,'-b')
plot(Phi_s,Yt4,'-r')




figure(2); clf; hold on;
plot(Phi_s,Yt1,'LineWidth',2)
plot(Phi_s,Yt2,'LineWidth',2)
plot(Phi_s,Yt3,'LineWidth',2)
plot(Phi_s,Yt4,'LineWidth',2)
legend('support=5','support=10','support=15','support=20')
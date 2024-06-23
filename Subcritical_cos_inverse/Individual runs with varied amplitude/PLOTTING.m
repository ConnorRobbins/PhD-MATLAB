data1=load('amp0.02.mat');
data2=load('amp0.05.mat');
data3=load('amp0.07.mat');
data4=load('amp0.1.mat');

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
legend('amp=0.02','amp=0.05','amp=0.07','amp=0.1')
ylabel('Yt')
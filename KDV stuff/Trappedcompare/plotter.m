load('fullynonlineartrapped.mat')
load('kdvtrapped.mat')

figure(1); clf; hold on;
plot(Phi,Ys_Newton,'-b')
plot(Phi,Yt_Newton,'-r')
plot(x,1+eta,'--k')
plot(x,forcing,'--k')
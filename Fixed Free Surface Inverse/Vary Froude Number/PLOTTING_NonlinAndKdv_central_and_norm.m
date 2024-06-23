clear

%%%%%% load the positive surface amplitude data
s=load('a1E-1_b3E-1_N641.mat');
[p1_F,p1_nonlin_norm,p1_nonlin_central,p1_kdv_norm,p1_kdv_central]=deal(s.Froudes,s.Yt_norm_store,s.Yt_0_store,s.kdv_norm_store,s.kdv_0_store);
s=load('a1p5E-1_b3E-1_N641.mat');
[p2_F,p2_nonlin_norm,p2_nonlin_central,p2_kdv_norm,p2_kdv_central]=deal(s.Froudes,s.Yt_norm_store,s.Yt_0_store,s.kdv_norm_store,s.kdv_0_store);
s=load('a2E-1_b3E-1_N641.mat');
[p3_F,p3_nonlin_norm,p3_nonlin_central,p3_kdv_norm,p3_kdv_central]=deal(s.Froudes,s.Yt_norm_store,s.Yt_0_store,s.kdv_norm_store,s.kdv_0_store);

%%%%%%% load the negative surface amplitude data
s=load('a-1E-1_b3E-1_N641.mat');
[n1_F,n1_nonlin_norm,n1_nonlin_central,n1_kdv_norm,n1_kdv_central]=deal(s.Froudes,s.Yt_norm_store,s.Yt_0_store,s.kdv_norm_store,s.kdv_0_store);
s=load('a-1p5E-1_b3E-1_N641.mat');
[n2_F,n2_nonlin_norm,n2_nonlin_central,n2_kdv_norm,n2_kdv_central]=deal(s.Froudes,s.Yt_norm_store,s.Yt_0_store,s.kdv_norm_store,s.kdv_0_store);
s=load('a-2E-1_b3E-1_N641.mat');
[n3_F,n3_nonlin_norm,n3_nonlin_central,n3_kdv_norm,n3_kdv_central]=deal(s.Froudes,s.Yt_norm_store,s.Yt_0_store,s.kdv_norm_store,s.kdv_0_store);




orange=[0.9290 0.6940 0.1250];
blue=[0 0.4470 0.7410];
red=[0.8500 0.3250 0.0980];
black=[0 0 0];

%%%%%% PLOT  positive surf norms
figure(1); clf; hold on; box on;
plot(p1_F,p1_nonlin_norm,DisplayName='$\alpha=0.1$',LineWidth=2,Color=red)
plot(p2_F,p2_nonlin_norm,DisplayName='$\alpha=0.15$',LineWidth=2,Color=blue)
plot(p3_F,p3_nonlin_norm,DisplayName='$\alpha=0.2$',LineWidth=2,Color=black)
%
plot(p1_F,p1_kdv_norm,LineWidth=1,Color=red,HandleVisibility="off",LineStyle="--")
plot(p2_F,p2_kdv_norm,LineWidth=1,Color=blue,HandleVisibility="off",LineStyle="--")
plot(p3_F,p3_kdv_norm,LineWidth=1,Color=black,HandleVisibility="off",LineStyle="--")
legend('Interpreter','Latex','FontSize',14)
ylabel('$\|y_b\|$',Interpreter='latex',FontSize=18)
xlabel('$F$',Interpreter='latex',FontSize=18)
ylim([0,0.4])



%%%%%% PLOT  negative surf norms
figure(2); clf; hold on; box on;
plot(n1_F,n1_nonlin_norm,DisplayName='$\alpha=-0.1$',LineWidth=2,Color=red)
plot(n2_F,n2_nonlin_norm,DisplayName='$\alpha=-0.15$',LineWidth=2,Color=blue)
plot(n3_F,n3_nonlin_norm,DisplayName='$\alpha=-0.2$',LineWidth=2,Color=black)
%
plot(n1_F,n1_kdv_norm,LineWidth=1,Color=red,HandleVisibility="off",LineStyle="--")
plot(n2_F,n2_kdv_norm,LineWidth=1,Color=blue,HandleVisibility="off",LineStyle="--")
plot(n3_F,n3_kdv_norm,LineWidth=1,Color=black,HandleVisibility="off",LineStyle="--")
%
legend('Interpreter','Latex','FontSize',14,Location='northwest')
ylabel('$\|y_b\|$',Interpreter='latex',FontSize=18)
xlabel('$F$',Interpreter='latex',FontSize=18)
ylim([0,0.25])


%%%%%% PLOT  positive surf centralvals
figure(3); clf; hold on; box on;
plot(p1_F,p1_nonlin_central,DisplayName='$\alpha=0.1$',LineWidth=2,Color=red)
plot(p2_F,p2_nonlin_central,DisplayName='$\alpha=0.15$',LineWidth=2,Color=blue)
plot(p3_F,p3_nonlin_central,DisplayName='$\alpha=0.2$',LineWidth=2,Color=black)
%
plot(p1_F,p1_kdv_central,LineWidth=1,Color=red,HandleVisibility="off",LineStyle="--")
plot(p2_F,p2_kdv_central,LineWidth=1,Color=blue,HandleVisibility="off",LineStyle="--")
plot(p3_F,p3_kdv_central,LineWidth=1,Color=black,HandleVisibility="off",LineStyle="--")
%
legend('Interpreter','Latex','FontSize',14,Location='northwest')
ylabel('$y_b(0)$',Interpreter='latex',FontSize=18)
xlabel('$F$',Interpreter='latex',FontSize=18)
ylim([-0.3,0.1])




%%%%%% PLOT  negative surf centralvals
figure(4); clf; hold on; box on;
plot(n1_F,n1_nonlin_central,DisplayName='$\alpha=-0.1$',LineWidth=2,Color=red)
plot(n2_F,n2_nonlin_central,DisplayName='$\alpha=-0.15$',LineWidth=2,Color=blue)
plot(n3_F,n3_nonlin_central,DisplayName='$\alpha=-0.2$',LineWidth=2,Color=black)
%
plot(n1_F,n1_kdv_central,LineWidth=1,Color=red,HandleVisibility="off",LineStyle="--")
plot(n2_F,n2_kdv_central,LineWidth=1,Color=blue,HandleVisibility="off",LineStyle="--")
plot(n3_F,n3_kdv_central,LineWidth=1,Color=black,HandleVisibility="off",LineStyle="--")
%
legend('Interpreter','Latex','FontSize',14,Location='northeast')
ylabel('$y_b(0)$',Interpreter='latex',FontSize=18)
xlabel('$F$',Interpreter='latex',FontSize=18)
ylim([-0.15,0.05])
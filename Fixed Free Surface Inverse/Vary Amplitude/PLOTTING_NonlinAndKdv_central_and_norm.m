clear

%%%%%% load the supercritical data
s=load('F1p1_loop_NegAndPos.mat');
[p1_amps,p1_nonlin_norm,p1_nonlin_central,p1_kdv_norm,p1_kdv_central]=deal(s.amplitudes,s.Yt_norm_store,s.Yt_0_store,s.kdv_norm_store,s.kdv_0_store);


%%%%%%% load the subcritical data
%s=load('F0p9_loop_NegAndPos.mat'); % set ylims to ([-0.25,0.02])
s=load('F0p97_loop_NegAndPos.mat');  % set ylims to ([-0.25,0.02])
[n1_amps,n1_nonlin_norm,n1_nonlin_central,n1_kdv_norm,n1_kdv_central]=deal(s.amplitudes,s.Yt_norm_store,s.Yt_0_store,s.kdv_norm_store,s.kdv_0_store);



orange=[0.9290 0.6940 0.1250];
blue=[0 0.4470 0.7410];
red=[0.8500 0.3250 0.0980];
black=[0 0 0];


%%%%%% PLOT  supercritical norms
figure(1); clf; hold on; box on;
plot(p1_amps,p1_nonlin_norm,DisplayName='$\alpha=0.1$',LineWidth=2,Color=black)
%
plot(p1_amps,p1_kdv_norm,LineWidth=1,Color=black,HandleVisibility="off",LineStyle="--")
%legend('Interpreter','Latex','FontSize',14)
ylabel('$\|y_b\|$',Interpreter='latex',FontSize=18)
xlabel('$\alpha$',Interpreter='latex',FontSize=18)
%ylim([0,0.4])




%%%%%% PLOT  supercritical centralvals
figure(2); clf; hold on; box on; grid on;
plot(p1_amps,p1_nonlin_central,DisplayName='$\alpha=0.1$',LineWidth=2,Color=black)
%
plot(p1_amps,p1_kdv_central,LineWidth=1,Color=black,HandleVisibility="off",LineStyle="--")
%
%legend('Interpreter','Latex','FontSize',14,Location='northwest')
ylabel('$y_b(0)$',Interpreter='latex',FontSize=18)
xlabel('$\alpha$',Interpreter='latex',FontSize=18)
%ylim([-0.3,0.1])






%%%%%% PLOT  subcritical norms
figure(3); clf; hold on; box on;
plot(n1_amps,n1_nonlin_norm,DisplayName='$\alpha=-0.1$',LineWidth=2,Color=black)
%
plot(n1_amps,n1_kdv_norm,LineWidth=1,Color=black,HandleVisibility="off",LineStyle="--")
%
%legend('Interpreter','Latex','FontSize',14,Location='northwest')
ylabel('$\|y_b\|$',Interpreter='latex',FontSize=18)
xlabel('$\alpha$',Interpreter='latex',FontSize=18)
%ylim([0,0.25])






%%%%%% PLOT  subcritical centralvals
figure(4); clf; hold on; box on; grid on;
plot(n1_amps,n1_nonlin_central,DisplayName='$\alpha=-0.1$',LineWidth=2,Color=black)
%
plot(n1_amps,n1_kdv_central,LineWidth=1,Color=black,HandleVisibility="off",LineStyle="--")
%
%legend('Interpreter','Latex','FontSize',14,Location='northeast')
ylabel('$y_b(0)$',Interpreter='latex',FontSize=18)
xlabel('$\alpha$',Interpreter='latex',FontSize=18)
%ylim([-0.15,0.05])
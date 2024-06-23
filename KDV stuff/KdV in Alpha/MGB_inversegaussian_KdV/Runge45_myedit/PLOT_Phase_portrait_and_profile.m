s1=load('alpha1p25result.mat');
s2=load('alpha1p5result.mat');
s3=load('alpha3result.mat');
s4=load('alpha5result.mat');



algebraic_x=linspace(-1.1,0,100);
algebraic_y=sqrt((-2/3)*algebraic_x.^3);
forcing_x=linspace(-10,10,1000);


figure(1); clf; hold on; box on;
plot(s1.u(1,:),s1.u(2,:),Color=[150,0,150]./255,LineWidth=2,DisplayName="$\alpha=1.25$")
plot(s2.u(1,:),s2.u(2,:),Color=[210,160,80]./255,LineWidth=2,DisplayName="$\alpha=1.5$")
plot(s3.u(1,:),s3.u(2,:),Color=[62,150,81]./255,LineWidth=2,DisplayName="$\alpha=3$")
plot(s4.u(1,:),s4.u(2,:),Color=[40,40,40]./255,LineWidth=2,DisplayName="$\alpha=5$")
plot(algebraic_x,algebraic_y,'--r',LineWidth=2,HandleVisibility='off')
plot(algebraic_x,-algebraic_y,'--r',LineWidth=2,HandleVisibility='off')
legend(Interpreter="latex",FontSize=14,Location="southwest")
ylabel('$\frac{\mathrm{d}u}{\mathrm{d}\xi}$',Interpreter='latex',FontSize=20,Rotation=0)
xlabel('$u$',Interpreter='latex',FontSize=18)
xlim([-1.6,2.1])
ylim([-3.5,1])


figure(2); clf; hold on; box on;
plot(s1.xx,s1.uu,Color=[150,0,150]./255,LineWidth=2,DisplayName="$\alpha=1.25$")
plot(-s1.xx,s1.uu,Color=[150,0,150]./255,LineWidth=2,HandleVisibility='off')
plot(s2.xx,s2.uu,Color=[210,160,80]./255,LineWidth=2,DisplayName="$\alpha=1.5$")
plot(-s2.xx,s2.uu,Color=[210,160,80]./255,LineWidth=2,HandleVisibility='off')
plot(s3.xx,s3.uu,Color=[62,150,81]./255,LineWidth=2,DisplayName="$\alpha=3$")
plot(-s3.xx,s3.uu,Color=[62,150,81]./255,LineWidth=2,HandleVisibility='off')
plot(s4.xx,s4.uu,Color=[40,40,40]./255,LineWidth=2,DisplayName="$\alpha=5$")
plot(-s4.xx,s4.uu,Color=[40,40,40]./255,LineWidth=2,HandleVisibility='off')
legend(Interpreter="latex",FontSize=14,Location="northeast")
xlim([-10,10])
ylabel('$u(\xi)$',Interpreter='latex',FontSize=18)
xlabel('$\xi$',Interpreter='latex',FontSize=18)


figure(3); clf; hold on; box on;
plot(forcing_x,forcing_for_plot(forcing_x,1.25),Color=[150,0,150]./255,LineWidth=2,DisplayName="$\alpha=1.25$")
plot(forcing_x,forcing_for_plot(forcing_x,1.5),Color=[210,160,80]./255,LineWidth=2,DisplayName="$\alpha=1.5$")
plot(forcing_x,forcing_for_plot(forcing_x,3),Color=[62,150,81]./255,LineWidth=2,DisplayName="$\alpha=3$")
plot(forcing_x,forcing_for_plot(forcing_x,5),Color=[40,40,40]./255,LineWidth=2,DisplayName="$\alpha=5$")
legend(Interpreter="latex",FontSize=14,Location="northeast")
% xlim([-10,10])
ylabel('$\alpha s(\xi)$',Interpreter='latex',FontSize=18)
xlabel('$\xi$',Interpreter='latex',FontSize=18)



function forcing_y=forcing_for_plot(forcing_x,alpha)
    forcing_y=alpha*((((4*forcing_x.^2)-2).*exp(-(forcing_x.^2))) + exp(-2*(forcing_x.^2)));
end
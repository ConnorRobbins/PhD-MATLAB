s1=load('forced_solitary_a2E-2_F1p1_N641.mat');
s2=load('forced_uniformstream_a2E-2_F1p1.mat');
s3=load('unforced_solitary_wave_F1p1.mat');

Phi=s3.Phi;


% Unforced case
figure(1); clf; hold on; box on;
plot(Phi+0.75,s3.Ys_Newton,'-k','LineWidth',2)
plot(Phi,ones(641,1),'LineWidth',2)
xlim([-20,20])
ylim([0.95,1.25])
ylabel('$y_f(\phi)$',Interpreter='latex',FontSize=18)
xlabel('$\phi$',Interpreter='latex',FontSize=18)

% Forced case
figure(2);clf; hold on; box on;
plot(Phi,s1.Ys_Newton,'-k','LineWidth',2)
plot(Phi,s2.Ys_Newton,'LineWidth',2)
ylim([0.95,1.25])
ylabel('$y_f(\phi)$',Interpreter='latex',FontSize=18)
xlabel('$\phi$',Interpreter='latex',FontSize=18)



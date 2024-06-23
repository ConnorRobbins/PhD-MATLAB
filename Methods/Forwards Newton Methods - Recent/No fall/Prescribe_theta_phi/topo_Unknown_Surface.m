

L_width = 100; %Phi range on either side of the origin
N_points = 341;
Froude = 0.7;

ijac = 1; %`1= calculate jac every newton loop, 0=only first loop

Phi = linspace(-L_width,L_width,N_points);

%% Forcings
Pressure_amplitude = -0.00;
Topo_amplitude=0.01; %% based on a*exp(-(bx)^2)
Topo_b_width=0.3; %% based on a*exp(-(bx)^2)
P = Pressure_amplitude*exp(-(Phi).^2); % pressure distribution
Theta_bottom=atan(-2*Topo_amplitude*(Topo_b_width^2)*Phi.*exp(-(Topo_b_width*Phi).^2)); %theta on topography

%% Initial guess for Theta
Theta_initial_guess=zeros(N_points,1);


%% Run Newton function

[Phi,Xs,Ys,Xt,Yt,Theta_surface,Theta_bottom]=surfacefinder_topography_fixedFP(L_width,N_points,Froude,ijac,P,Theta_bottom,Theta_initial_guess);

%%

figure
subplot(2,1,1)
plot(Xs,Ys-1,'-ob')
hold on
plot(Xt,Yt,'-or')
legend('Surface-1','Topography')
xlabel('X')
ylabel('y')
xlim([Xt(1) Xt(end)])

subplot(2,1,2)
plot(Xs,Theta_surface,'-xb')
hold on
plot(Xt,Theta_bottom,'-xr')
legend('Surface','Topography')
xlabel('X')
ylabel('theta')
xlim([Xt(1) Xt(end)])





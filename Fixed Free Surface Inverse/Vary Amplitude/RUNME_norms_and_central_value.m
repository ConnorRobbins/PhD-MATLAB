clear
saveplots=0;

%% Physical parameters 


N=641; L=20; b_width=0.3; Froude=0.97; %amplitudes=-[0.05,0.1,0.125,0.15,0.2]; %truncation_ranks=[];


%amplitudes=linspace(0.001,0.2,500);
amplitudes=linspace(-0.01,0.01,501);

truncation_ranks=[54];


[N_t,N_s]=deal(N);
Phi_s=linspace(-L,L,N_s); Phi_t=linspace(-L,L,N_t);
P=zeros(1,N_s);








[N_t,N_s]=deal(N);
Phi_s=linspace(-L,L,N_s); Phi_t=linspace(-L,L,N_t);
P=zeros(1,N_s);



loopnumber=numel(amplitudes);
steps=num2str(loopnumber);


Yt_norm_store=zeros(1,loopnumber);
P_norm_store=zeros(1,loopnumber);
Yt_0_store=zeros(1,loopnumber);
P_0_store=zeros(1,loopnumber);
%kdv_0_store=zeros(1,loopnumber);


if numel(truncation_ranks)==1
    truncation_ranks=ones(1,loopnumber)*truncation_ranks;
end

for i=1:loopnumber
    %
    disp("Starting step "+num2str(i)+"/"+steps)
    amplitude=amplitudes(i);
    truncation_rank=truncation_ranks(i);
    Ys=1+amplitude*exp(-(b_width*Phi_s).^2);
    Theta_approx_for_pressure=atan(-2*amplitude*(b_width^2)*Phi_s.*exp(-(b_width*Phi_s).^2))';

    %% %% Perform the topography calculations
    %Form the truncated pseudoinverse
    [truncated_psuedoinv] = FUN_make_SVD_matrix(L,N_s,N_t,Froude,truncation_rank);
    % form the RHS of the matrix eqn and output Theta_f
    [D,Theta] = FUN_make_SVD_RHS(Ys,P,Froude,L,N_s,N_t);
    %Solve the truncated system
    Theta_b=truncated_psuedoinv*D;
    [Yt,~] = Y_after_SVD_fun(L,N_t,N_s,Theta_b,Theta);
    Yt_norm_store(i)=sqrt(trapz(Phi_t,Yt.^2));
    Yt_0_store(i)=Yt((N_t+1)/2);


    %% %% Perform the pressure calculations
    [~,~,~,~,~,~,~,inverse_Pressure] = FUN_Inverse_Pressure_NO_TOPOGRAPHY(L,N_s,N_t,Froude,Theta_approx_for_pressure);
    P_norm_store(i)=sqrt(trapz(Phi_s,inverse_Pressure.^2));
    P_0_store(i)=inverse_Pressure((N_s+1)/2);


    %% %% Calculate the central value of the kdv forcing
    eta=Ys-1;
    eta_xx=[0,(eta(3:end)-2*eta(2:end-1)+eta(1:end-2))/((Phi_s(2)-Phi_s(1)).^2),0];
    kdv_forcing=2*(Froude-1)*eta -(3*(eta.^2)/2) -(eta_xx/3);
    kdv_0_store(i)=kdv_forcing((N_s+1)/2);
    kdv_trapz_norm_store(i)=sqrt(trapz(Phi_s,kdv_forcing.^2));


end

kdv_norm_store= abs(amplitudes*(pi^(1/4))).*sqrt(((4*(Froude-1).^2)/(b_width*sqrt(2)))    + ((4*b_width*(Froude-1) + b_width^3)/(3*sqrt(2)))  -  ((6*amplitudes*(Froude-1))/(b_width*sqrt(3)))   -   (4*amplitudes*b_width/(3*sqrt(3))) + ((9*amplitudes.^2)/(8*b_width)));
kdv_0_store=amplitudes.*(2*(Froude-1) + (2*(b_width^2)/3)-(3*amplitudes/2));



%%





figure(1); clf; hold on; box on;
plot(amplitudes,P_norm_store,':k','LineWidth',1.5,DisplayName='$\|P\|$: Fully Nonlinear')
plot(amplitudes,Yt_norm_store,'--k',DisplayName='$\|y_b\|$: Fully Nonlinear')
plot(amplitudes,kdv_norm_store,'-k',DisplayName='$\|f\|$: KdV')
legend('Interpreter','Latex','FontSize',11)
xlabel('$\alpha$','Interpreter','Latex','FontSize',16)
ylabel('$\|$Forcing$\|$','Interpreter','Latex','FontSize',16)
% xlim([0.95,1.25])
% ylim([0,0.2])
% axes('Position',[.3 .7 .2 .2]); box on; hold on;
% plot(amplitudes,P_norm_store,':k','LineWidth',1.5)
% plot(amplitudes,Yt_norm_store,'--k')
% plot(amplitudes,kdv_norm_store,'-k')
% xlim([1.001,1.006])
% ylim([0.084,0.088])




figure(2);  clf; hold on; box on;
plot(amplitudes,P_0_store,':k','LineWidth',2)
plot(amplitudes,Yt_0_store,'--k')
plot(amplitudes,kdv_0_store,'-k')
xlabel('$\alpha$','Interpreter','Latex','FontSize',16)
ylabel('Forcing central value','FontSize',16,'Interpreter','Latex')
legend('$P(0)$: Fully Nonlinear',' $y_b(0)$: Fully Nonlinear','$f(0)$: KdV','Interpreter','Latex','FontSize',11)
% ylim([-0.4,0.2])
% xlim([0.8,1.3])
% axes('Position',[.3 .7 .2 .2]); box on; hold on;
% plot(amplitudes,P_0_store,':k','LineWidth',2)
% plot(amplitudes,Yt_0_store,'--k')
% plot(amplitudes,kdv_0_store,'-k')
% xlim([1.0075,1.025])
% ylim([-0.048,-0.037])

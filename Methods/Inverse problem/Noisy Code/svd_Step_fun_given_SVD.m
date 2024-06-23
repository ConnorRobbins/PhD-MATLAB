function [Phi_t,Phi_s,Theta,Theta_b,D] = svd_Step_fun_given_SVD(L,Froude,N_s,N_t,P,Ys,truncated_inv)
%SVD_NO_P_FUN Summary of this function goes here
%   Detailed explanation goes here


% Make Ys a row
Ys=reshape(Ys,1,N_s);


% Create meshes
Phi_t=linspace(-L,L,N_t);
Phi_s=linspace(-L,L,N_s);
%del_Phi_t=2*L/(N_t-1);
del_Phi_s=2*L/(N_s-1);
%Phi_t_mid=0.5*(Phi_t(1:end-1)+Phi_t(2:end));
Phi_s_mid=0.5*(Phi_s(1:end-1)+Phi_s(2:end));




%% Calculate Theta from Ys

%Calculate dy/dphi by centered difference
dydphi=zeros(1,N_s);
dydphi(2:end-1)=0.5*(Ys(3:end)-Ys(1:end-2))/del_Phi_s;

%Calculate tau_f from Bernoulli's eqn
tau_f= 0.5*log(1+(2*(1-Ys-P)/(Froude^2)));


Theta=asin(dydphi.*exp(tau_f));


%% Form RHS of matrix equation

% Initialise RHS vector
D=zeros(N_s,1);

% Form surface variables
Ys_mid=0.5*(Ys(1:end-1)+Ys(2:end));
P_mid=0.5*(P(1:end-1)+P(2:end));


% D=surftrapz-logterm
G_k_I=Theta'./(1-exp(pi*(Phi_s_mid-Phi_s')));
surface_trapz= del_Phi_s*((0.5* (G_k_I(1,:) + G_k_I(N_s,:) ) ) + sum(G_k_I(2:N_s-1,:)))';
log_term= 0.5*log(1 + ( (2/(Froude^2))*(1-Ys_mid - P_mid) )   )';

D(1:end-1)=surface_trapz-log_term;





%% Perform the multiplication



Theta_b=truncated_inv*D;


end


function [Phi_t,Phi_s,Theta,Theta_b_coeff_matrix,Theta_b,U,singularvals,V,D,maxrank] = svd_Step_fun(L,Froude,N_s,N_t,P,Ys)
%svd_Step_fun : Solving the equation Theta_b_coeff_matrix * theta_b = D
%Computes the SVD of the coefficient matrix then uses this to
%calculate theta_b at every possible truncation rank


% Ensure Ys is a row vector
Ys=reshape(Ys,1,N_s);
P=reshape(P,1,N_s);


% Create meshes
Phi_t=linspace(-L,L,N_t);
Phi_s=linspace(-L,L,N_s);
del_Phi_t=2*L/(N_t-1);
del_Phi_s=2*L/(N_s-1);
Phi_t_mid=0.5*(Phi_t(1:end-1)+Phi_t(2:end));
Phi_s_mid=0.5*(Phi_s(1:end-1)+Phi_s(2:end));



%% Form coefficient matrix
Theta_b_coeff_matrix=zeros(N_s,N_t);

% Form the matrix of 1/1+exp() 
middle_mat=Phi_s_mid'-Phi_t;
middle_mat=1+exp(pi*middle_mat);
middle_mat=1./middle_mat;

% Fill first and last columns of coeff matrix from exponentials
Theta_b_coeff_matrix(1:N_s-1,1)=0.5*del_Phi_t*middle_mat(:,1);
Theta_b_coeff_matrix(1:N_s-1,N_t)=0.5*del_Phi_t*middle_mat(:,N_t);
% Fill middle columns
Theta_b_coeff_matrix(1:N_s-1,2:N_t-1)=del_Phi_t*middle_mat(:,2:N_t-1);


% Fill final row of coeff matrix with BC
if Froude<1
    %Theta_b_coeff_matrix(N_s,N_t)=1;  % The 0 is the value of theta at last grid point (BC for subcritical)
    Theta_b_coeff_matrix(N_s,1)=1; 
else
    Theta_b_coeff_matrix(N_s,1)=1;  % The 0 is the value of theta at first grid point (BC for supercritical)
end





%% Calculating Theta (on surface) from Ys 

%Calculate dy/dphi by centered difference, setting first and last value to
%0
dydphi=zeros(1,N_s);
dydphi(2:end-1)=0.5*(Ys(3:end)-Ys(1:end-2))/del_Phi_s;

%Calculate tau_f from rearranged Bernoulli's eqn
tau_f= 0.5*log(1+(2*(1-Ys-P)/(Froude^2)));

%Now evaluate Theta
Theta=asin(dydphi.*exp(tau_f));


%% Form RHS of matrix equation

% Initialise RHS vector
D=zeros(N_s,1);

% Form surface variables
Ys_mid=0.5*(Ys(1:end-1)+Ys(2:end));
P_mid=0.5*(P(1:end-1)+P(2:end));

%surftrapz is the integral over the surface and logterm is the log part of
%the RHS (D) of matrix equation
% D=surftrapz-logterm
G_k_I=Theta'./(1-exp(pi*(Phi_s_mid-Phi_s')));
surface_trapz= del_Phi_s*((0.5* (G_k_I(1,:) + G_k_I(N_s,:) ) ) + sum(G_k_I(2:N_s-1,:)))';
log_term= 0.5*log(1 + ( (2/(Froude^2))*(1-Ys_mid - P_mid) )   )';

D(1:end-1)=surface_trapz-log_term; %final entry is a 0 to correspond to BC row of matrix



%% Perform the SVD and truncate at every rank to calculate all possible truncated solutions


maxrank=min(N_s,N_t); %Identify maxrank
Theta_b=zeros(N_t,maxrank); %Initialise Theta_b storage, the kth column will be the solution at rank k

%[U,singularvals,V]=svd(Theta_b_coeff_matrix,'vector'); %Calculate SVD
[U,singularvals,V]=svd(Theta_b_coeff_matrix); %Calculate SVD
singularvals=diag(singularvals);

u=U';

truncated_inv=zeros(N_t,N_s);
% Now loop over the ranks 
for k=1:maxrank
    truncated=(1/singularvals(k))*V(:,k)*u(k,:);
    truncated_inv=truncated_inv+truncated;
    Theta_b(:,k)=truncated_inv*D;
end



end


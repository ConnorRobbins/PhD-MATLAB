function [truncated_inv] = FUN_make_SVD_matrix(L,N_s,N_t,Froude,truncation_rank)
%FUN_MAKE_SVD_MATRIX Summary of this function goes here
%   Detailed explanation goes here


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


% % Fill final row of coeff matrix with BC
% if Froude<1
%     Theta_b_coeff_matrix(N_s,N_t)=1;  % The 0 is the value of theta at last grid point (BC for subcritical)
%     %Theta_b_coeff_matrix(N_s,1)=1; 
% else
%     Theta_b_coeff_matrix(N_s,1)=1;  % The 0 is the value of theta at first grid point (BC for supercritical)
% end
Theta_b_coeff_matrix(N_s,1)=1;  % The 0 is the value of theta at first grid point (BC for supercritical)





%% Now truncate the matrix



[U,singularvals,V]=svd(Theta_b_coeff_matrix,'vector'); %Calculate SVD
% [U,singularvals,V]=svd(Theta_b_coeff_matrix); %Calculate SVD
% singularvals=diag(singularvals);

u=U';

truncated_inv=zeros(N_t,N_s);
% Now loop over the ranks 
for k=1:truncation_rank
    truncated=(1/singularvals(k))*V(:,k)*u(k,:);
    truncated_inv=truncated_inv+truncated;
end












end


function [Phi_t,Phi_s,Theta,Theta_b_coeff_matrix,Theta_b,U,singularvals,V,D,maxrank,Theta_b_norm,rankSVD,Yt,Ys_after]=func_unperturbed_problem(L,Froude,N_s,N_t,P,Ys,boo_auto_pick_rank,rankSVD,normmax,grad_tol)

% Call svd function to calculate SVD and compute theta_b at all truncation
% ranks
[Phi_t,Phi_s,Theta,Theta_b_coeff_matrix,Theta_b,U,singularvals,V,D,maxrank] = svd_Step_fun(L,Froude,N_s,N_t,P,Ys);

% Calculate the norm of each rank's solution for theta_b, this will be used
% to help decide at which rank we truncate
Theta_b_norm=zeros(1,maxrank);
for i=1:maxrank
    Theta_b_norm(i)=sqrt(trapz(Phi_t,Theta_b(:,i).^2));
end

% If rank selection for unperturbed problem is on then call funciton to
% find rank
if boo_auto_pick_rank==1
    [rankSVD] = gradient_rank_selector_capped(Theta_b_norm,grad_tol,normmax) % Select rank for TSVD based on gradient of norm (has built in: will not allow norm of theta_b at selected to be above 10^normmax)
end

[Yt,Ys_after] = Y_after_SVD_fun(L,N_t,N_s,Theta_b(:,rankSVD),Theta); %Calculate Yt from theta_b at selected truncation rank

%%%%% This solves the inverse problem by TSVD once without perturbation and
%%%%% then (adds random noise and solves it again) "Random_error_Number"
%%%%% times. The specified number of end points do not have noise added to them.
%%

clear


load('fall_soln_a0.1_N721.mat')
N_s=N; N_t=N;

%% Error based parameters

grad_tol=1E-3;  %Tolerance used to decide which parts of solution_norm-rank curve are "flat enough"
normmax=4; % If the norm of the theta_b is above 10^normmax at rank k then disallow the selection of k as the truncation rank


boo_auto_pick_rank=0; %If this is set to 1 the code will ignore rankSVD set below and instead decide on rankSVD based on going 80% along the flat part of the norm curve 
rankSVD=90; %WILL BE OVERWRITTEN IF auto_pick_unperturbed_rank=1;   The fixed truncation rank for applying TSVD to the unperturbed problem






%% Solve the unperturbed problem

[Phi_t,Phi_s,Theta,Theta_b_coeff_matrix,Theta_b,U,singularvals,V,D,maxrank,Theta_b_norm,rankSVD,Yt,Ys_after]=func_unperturbed_problem(L,Froude,N_s,N_t,P,Ys_Newton,boo_auto_pick_rank,rankSVD,normmax,grad_tol);






%% Plots 


% figure(4); clf; hold on;
% plot(Phi_t,Theta_b_store)
% ylabel('$\theta_b$','Interpreter','latex','Rotation',0,'FontSize',16)




figure(9); clf; hold on;
plot(Phi_s,Ys_Newton,'-b')
plot(Phi_t,Yt,'-k')
ylabel('$y$','Interpreter','latex','Rotation',0,'FontSize',16)




figure(10); clf; hold on;
plot(Phi_t,Yt,'-k')
plot(Phi_t,Yt_Newton,'--r')
ylabel('Yb')
legend('Inverse','Original Newton')









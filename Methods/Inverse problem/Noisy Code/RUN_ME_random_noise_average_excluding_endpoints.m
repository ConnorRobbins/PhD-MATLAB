%%%%% This solves the inverse problem by TSVD once without perturbation and
%%%%% then (adds random noise and solves it again) "Random_error_Number"
%%%%% times. The specified number of end points do not have noise added to them.
%%

clear

%% Physical parameters 
L=25; %domain truncation
Froude=1.2;
amplitude=0.1;  %surface parameter, see form of Ys
b_width=0.3;    %surface parameter, see form of Ys
topo_centre=0;  %surface parameter, see form of Ys

N_t=151;         %Number of points on topography
N_s=N_t;        %Number of points on surface

Phi_t=linspace(-L,L,N_t); %Creates topography mesh
Phi_s=linspace(-L,L,N_s); %Creates surface mesh

P=zeros(1,N_s);     %Pressure on surface

Ys=1+amplitude*exp(-(b_width*(Phi_s-topo_centre)).^2);
%Ys=1+amplitude*cos(Phi_s).*exp(-(b_width*(Phi_s-topo_centre)).^2);

%% Error based parameters

grad_tol=1E-3;  %Tolerance used to decide which parts of solution_norm-rank curve are "flat enough"
Random_error_Number=100; %Number of times to perturb free surface and solve
normmax=4; % If the norm of the theta_b is above 10^normmax at rank k then disallow the selection of k as the truncation rank
Noise_range=L; %values between which noise is added to the surface

unperturbed_endpoints=2; %number of endpoints on each side without added noise (overrides noise support if needed)
random_noise_std=5E-3; % standard deviation of the normal dist. from which noise is drawn

auto_pick_unperturbed_rank=0; %If this is set to 1 the code will ignore rankSVD set below and instead decide on rankSVD based on going 80% along the flat part of the norm curve 

rankSVD=50; %WILL BE OVERWRITTEN IF auto_pick_unperturbed_rank=1;   The fixed truncation rank for applying TSVD to the unperturbed problem
rankSVD_noisy=40; %The fixed truncation rank for applying TSVD to the noisy problems





%% Solve the unperturbed problem

[Phi_t,Phi_s,Theta,Theta_b_coeff_matrix,Theta_b,U,singularvals,V,D,maxrank,Theta_b_norm,rankSVD,Yt,Ys_after]=func_unperturbed_problem(L,Froude,N_s,N_t,P,Ys,auto_pick_unperturbed_rank,rankSVD,normmax,grad_tol);



%%  Solve the perturbed problems


[Yt_std,Ys_std,Ys_noisy_store,Yt_noisy_store,Theta_b_store,truncated_inv,Yt_noisy_average]=func_perturbed_surface_problem(L,Froude,N_s,N_t,P,Ys,Phi_s,Phi_t,rankSVD_noisy,U,singularvals,V,unperturbed_endpoints,Noise_range,Random_error_Number,random_noise_std);


%% Plots 


% figure(4); clf; hold on;
% plot(Phi_t,Theta_b_store)
% ylabel('$\theta_b$','Interpreter','latex','Rotation',0,'FontSize',16)




figure(9); clf; hold on;
plot(Phi_s,Ys,'-b')
plot(Phi_s,Ys_noisy_store(:,1),'-r')
ylabel('$Y_ f$','Interpreter','latex','Rotation',0,'FontSize',16)
legend('True surface','Example of perturbed surface')



figure(11); clf; hold on;
plot(Phi_t,Yt,'-b')
plot(Phi_t,Yt_noisy_average,'-r')
legend('not perturbed','averaged')
ylabel('Yb')

%

std_no=2;
figure(12); clf; hold on;
%%% Plot unpeturbed solution
plot(Phi_t,Yt,'-k','LineWidth',1.2)
%%% Plot average solution
plot(Phi_t,Yt_noisy_average,'-','LineWidth',1.2,'Color',[0.5,0.5,0.5])
%%% plot noisy results std on top of average
plot(Phi_t,Yt_noisy_average+std_no*Yt_std,'--k','LineWidth',1)
plot(Phi_t,Yt_noisy_average-std_no*Yt_std,'--k','LineWidth',1)






%legend('Unperturbed Yt','Noisy avg Yt','Noisy avg Yt +- 2*std(topo)','','Unperturbed Yt +- 2*std(Ys)','')
legend('Unperturbed Yt','Noisy avg Yt','Noisy avg Yt +- 2*std(topo)','')


%ylim(ylimits)
xlabel('$\phi$','Interpreter','latex','FontSize',13)
ylabel('$y_b$','Interpreter','latex','FontSize',17)%,'Rotation',0)
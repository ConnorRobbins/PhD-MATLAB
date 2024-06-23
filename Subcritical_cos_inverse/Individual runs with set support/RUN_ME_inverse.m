%%%%% This solves the inverse problem by TSVD once without perturbation.
% The specified number of end points do not have noise added to them.
%%

clear

%% Physical parameters 
L=30; %domain truncation
Froude=0.8;
amplitude=0.01;  %surface parameter, see form of Ys
b_width=0.3;    %surface parameter, see form of Ys
topo_centre=15;  %surface parameter, see form of Ys


N_t=851;         %Number of points on topography
N_s=N_t;        %Number of points on surface

Phi_t=linspace(-L,L,N_t); %Creates topography mesh
Phi_s=linspace(-L,L,N_s); %Creates surface mesh

P=zeros(1,N_s);     %Pressure on surface

%Ys=1+amplitude*exp(-(b_width*(Phi_s-topo_centre)).^2);
%Ys=1+amplitude*cos(Phi_s).*exp(-(b_width*(Phi_s-topo_centre)).^2);
%Ys=1+amplitude*cos(Phi_s);
Ys=1+amplitude*0.5*(tanh(topo_centre+Phi_s)+tanh(topo_centre-Phi_s)).*cos(Phi_s);

figure(3); hold on;
plot(Phi_s,Ys)
%% Error based parameters

grad_tol=1E-3;  %Tolerance used to decide which parts of solution_norm-rank curve are "flat enough"
normmax=4; % If the norm of the theta_b is above 10^normmax at rank k then disallow the selection of k as the truncation rank


boo_auto_pick_rank=0; %If this is set to 1 the code will ignore rankSVD set below and instead decide on rankSVD based on going 80% along the flat part of the norm curve 
rankSVD=400; %WILL BE OVERWRITTEN IF auto_pick_unperturbed_rank=1;   The fixed truncation rank for applying TSVD to the unperturbed problem






%% Solve the unperturbed problem

[Phi_t,Phi_s,Theta,Theta_b_coeff_matrix,Theta_b,U,singularvals,V,D,maxrank,Theta_b_norm,rankSVD,Yt,Ys_after]=func_unperturbed_problem(L,Froude,N_s,N_t,P,Ys,boo_auto_pick_rank,rankSVD,normmax,grad_tol);






%% Plots 


% figure(4); clf; hold on;
% plot(Phi_t,Theta_b_store)
% ylabel('$\theta_b$','Interpreter','latex','Rotation',0,'FontSize',16)




figure(9); clf; hold on;
plot(Phi_s,Ys,'-b')
plot(Phi_t,Yt,'-k')
ylabel('$y$','Interpreter','latex','Rotation',0,'FontSize',16)




figure(10); clf; hold on;
plot(Phi_t,Yt,'-k')
ylabel('Yb')









%% Physical parameters 
clear

s=load('N1241trappedc3.mat');
L=s.L;
Froude=s.Froude;
P=s.P;
Ys=s.Ys_Newton;
Yt_orig=s.Yt_Newton;
N_t=s.N;
N_s=s.N;




%[Phi_t,Phi_s,Theta_b_coeff_matrix,singularvals,singularvals_larger,Theta_b,U,S,V,D,maxrank,Theta] = svd_Step_fun(L,Froude,N_s,N_t,P,amplitude,b_width);
[Phi_t,Phi_s,Theta,Theta_b_coeff_matrix,Theta_b,U,singularvals,V,D,maxrank] = svd_Step_fun(L,Froude,N_s,N_t,P,Ys);

%%

Theta_rank_k_norm=zeros(N_t,1);
Theta_rank_k_eqn_difference_norm=zeros(N_t,1);

for k=1:N_t
    Theta_rank_k=Theta_b(:,k);
    Theta_rank_k_norm(k)=sqrt(Theta_rank_k'*Theta_rank_k);
    Theta_rank_k_eqn_difference=Theta_b_coeff_matrix*Theta_rank_k - D;
    Theta_rank_k_eqn_difference_norm(k)=sqrt(Theta_rank_k_eqn_difference'*Theta_rank_k_eqn_difference);
end


[~,mindex]=min(log10(Theta_rank_k_eqn_difference_norm))

%%

example_rank=75;

[Yt,~] = Y_after_SVD_fun(L,N_t,N_s,Theta_b(:,example_rank),Theta);

figure(1); clf; hold on;
plot(log10(Theta_rank_k_norm))
plot(example_rank,log10(Theta_rank_k_norm(example_rank)),'-xr',MarkerSize=8)
ylabel('log10(ThetaK norm)')
xlabel('rank K')

figure(2); clf; hold on;
plot(Theta_rank_k_eqn_difference_norm)
plot(example_rank,Theta_rank_k_eqn_difference_norm(example_rank),'-xr',MarkerSize=8)
ylabel('Ax-b norm')

figure(3); clf; hold on;
plot(log10(Theta_rank_k_eqn_difference_norm))
plot(example_rank,log10(Theta_rank_k_eqn_difference_norm(example_rank)),'-xr',MarkerSize=8)
ylabel('log10(Ax-b norm)')

figure(4); clf; hold on;
plot(log10(Theta_rank_k_eqn_difference_norm),log10(Theta_rank_k_norm))
plot(log10(Theta_rank_k_eqn_difference_norm(example_rank)),log10(Theta_rank_k_norm(example_rank)),'-rx',MarkerSize=8)
xlabel('log10(A thetaK - b norm)')
ylabel('log10(thetaK norm)')




figure(5); clf; hold on;
plot(Phi_t,Yt)
ylabel('Yb')

figure(6); clf; hold on;
plot(Phi_t,Yt_orig,'-r')
plot(Phi_t,Yt,'-b')
ylabel('Yb')
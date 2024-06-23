clear


%load('661hydrofall_a02.mat')
load('661hydrofall_a005.mat')
N_s=N;
N_t=N;






[Phi_t,Phi_s,Theta_inverse,Theta_b_coeff_matrix,singularvals,singularvals_larger,Theta_b_matrix_inverse,U,S,V,D,maxrank] = svd_Step_fun(L,Froude,N_s,N_t,P,Ys_Newton);


%% 
rankSVD=81;

Theta_bottom_inverse=Theta_b_matrix_inverse(:,rankSVD);

[Yt_inverse,Ys_inverse,Xt_inverse,Xs_inverse] = variables_after_SVD(L,N_t,N_s,Theta_bottom_inverse,Theta_inverse);



%



%%%%%%%%%%%%%%%% Theta_b
figure(1); clf; hold on;
plot(Phi_t,Theta_bottom_inverse,'-r')
plot(Phi,Theta_bottom_Newton,'--b')
xlabel('Phi')
ylabel('Theta bottom')
legend('Inverse','Newton')

figure(2); clf; hold on;
plot(Xt_inverse,Theta_bottom_inverse,'-r')
plot(Xt_Newton,Theta_bottom_Newton,'--b')
xlabel('Xt')
ylabel('Theta bottom')
legend('Inverse','Newton')



%%%%%%%%%% Yb
figure(3); clf; hold on;
plot(Phi_t,Yt_inverse,'-r')
plot(Phi,Yt_Newton,'--b')
xlabel('Phi')
ylabel('Yb')
legend('Inverse','Newton')

figure(4); clf; hold on;
plot(Xt_inverse,Yt_inverse,'-r')
plot(Xt_Newton,Yt_Newton,'--b')
xlabel('Xt')
ylabel('Yb')
legend('Inverse','Newton')



%%%%%%%%% Ys

figure(5); clf; hold on;
plot(Xs_Newton,Ys_Newton,'-k')
plot(Phi_s,Ys_inverse,'-r')
plot(Phi,Ys_Newton,'--b')
ylabel('Yf')
legend('Xf from Newton','Phi from inverse','Phi from Newton')


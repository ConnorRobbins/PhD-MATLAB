clear
load('661hydrofall_a02.mat')
N_s=N;
N_t=N;




XPhi_line_est=8;

XPhigradleft=(Xs_Newton(XPhi_line_est)-Xs_Newton(1))/(Phi(XPhi_line_est)-Phi(1));
leftline=Xs_Newton(1)+(Phi+L)*XPhigradleft;

XPhigradright=(Xs_Newton(end)-Xs_Newton(end-XPhi_line_est+1))/(Phi(end)-Phi(end-XPhi_line_est+1));
rightline=Xs_Newton(end)+(Phi-L)*XPhigradright;




figure(7); clf; hold on; 
plot(Phi,Xs_Newton,'-k')
plot(Phi,leftline,'--r','LineWidth',2)
plot(Phi,rightline,'--b','LineWidth',2)
ylabel('Xs')
xlabel('Phi')
legend('Solution','Linear approx at left','Linear approx at right')


YsPhiGuess = interp1(Xs_Newton,Ys_Newton,Phi);


%%
% Inverse from true Ys
[Phi_t,Phi_s,Theta_inverse_true,Theta_b_coeff_matrix_true,singularvals,singularvals_larger,Theta_b_matrix_inverse_true,U,S,V,D,maxrank] = svd_Step_fun(L,Froude,N_s,N_t,P,Ys_Newton);
% Inverse from true Ys
[Phi_t,Phi_s,Theta_inverse_guess,Theta_b_coeff_matrix_guess,singularvals,singularvals_larger,Theta_b_matrix_inverse_guess,U,S,V,D,maxrank] = svd_Step_fun(L,Froude,N_s,N_t,P,YsPhiGuess);

%% 
rankSVD=94;

Theta_bottom_inverse_true=Theta_b_matrix_inverse_true(:,rankSVD);
[Yt_inverse_true,Ys_inverse_true,Xt_inverse_true,Xs_inverse_true] = variables_after_SVD(L,N_t,N_s,Theta_bottom_inverse_true,Theta_inverse_true);

Theta_bottom_inverse_guess=Theta_b_matrix_inverse_guess(:,rankSVD);
[Yt_inverse_guess,Ys_inverse_guess,Xt_inverse_guess,Xs_inverse_guess] = variables_after_SVD(L,N_t,N_s,Theta_bottom_inverse_guess,Theta_inverse_guess);


figure(8); clf; hold on; 
plot(Phi_t,Theta_bottom_Newton,'-k')
plot(Phi_t,Theta_bottom_inverse_true,'-b')
plot(Phi_t,Theta_bottom_inverse_guess,'-r')
ylabel('Theta bottom')
xlabel('Phi')
legend('True','Inverse from true','Inverse from guess')



figure(9); clf; hold on;
plot(Xt_Newton,Yt_Newton,'-k')
plot(Xt_inverse_true,Yt_inverse_true,'-b')
plot(Xt_inverse_guess,Yt_inverse_guess,'-r')
ylabel('Yt')
xlabel('Xt')
legend('True','Inverse from true','Inverse from guess')




figure(6); clf; hold on;
plot(Phi,Ys_Newton,'-r')
plot(Xs_Newton,Ys_Newton,'-k')
plot(Phi_s,Ys_inverse_guess,'--g')
plot(Xs_inverse_guess,Ys_inverse_guess,'--b')
legend('Ys Phi Newt','Ys Xs Newt','Ys Phi Guess','Ys Xs Guess')



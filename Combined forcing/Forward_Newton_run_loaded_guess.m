clear 

N=641;
%filename='subcrit_guess';
filename='inverse_trial_N641_F0p8';



s=load(filename);
L=s.L;
Froude=s.Froude;
Phi_old=s.Phi;
Yt=s.Yt_example;
P=s.P;





Phi = linspace(-L,L,N);
ijac=1;





%initial guesses
Theta_init_guess=interp1(Phi_old,s.Theta,Phi);
Theta_b_init_guess=interp1(Phi_old,s.Theta_b(:,s.example_rank_1),Phi);





[Phi,Xs_Newton,Ys_Newton,Xt_Newton,Yt_Newton,Theta_Newton,Theta_bottom_Newton,Froude] = Newton_forward_fun(L,N,Froude,ijac,P,Theta_init_guess,Theta_b_init_guess,Yt);


%%

figure(11); clf; hold on;
plot(Phi,Ys_Newton,'-b')
plot(Phi,Yt_Newton,'-k')
plot(Phi,P,'--r')

%figure(13); clf;hold on;
% plot(Xs_Newton,Ys_Newton,'--b')
% plot(Xt_Newton,Yt_Newton,'--k')




figure(12); clf; hold on;
plot(Phi,Theta_bottom_Newton)
ylabel('Theta bottom')
xlabel('Phi')


figure(13); clf; hold on;
plot(Phi,1-Ys_Newton)
ylabel('eta')

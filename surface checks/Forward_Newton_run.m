clear 

%s=load("for_forward_test.mat");
s=load("subcheck_right.mat");

L=s.L;
N=s.N;


Phi=s.Phi_t;

%parameters
amplitude=s.amplitude;

F=s.Froude;
ijac=1;
P=s.P;



%initial guesses

Theta_init_guess=s.Theta;
Theta_b_init_guess=s.Theta_b(:,s.rankSVD);

Ys_orignal=s.Ys;



%[Phi,Xs_Newton,Ys_Newton,Xt_Newton,Yt_Newton,Theta_Newton,Theta_bottom_Newton,Froude] = Newton_forward_fun(L,N,Froude,amplitude,b_width,ijac,P,Theta_init_guess,Theta_b_init_guess);
[Phi,Xs,Ys,Xt,Yt,Theta,Theta_bottom,P] = surfacefinder_topography_fixedFP(L,N,F,ijac,P,Theta_b_init_guess,Theta_init_guess);
    

%%

figure(11); clf; hold on;
plot(Phi,Ys,'-b')
plot(Phi,Yt,'-k')
xlabel('Phi')

%figure(13); clf;hold on;
plot(Xs,Ys,'--b')
plot(Xt,Yt,'--k')
xlabel('x')




figure(12); clf; hold on;
plot(Phi,Theta_bottom)
ylabel('Theta bottom')
xlabel('Phi')


figure(14); clf; hold on;
plot(Phi,Ys,'-r')
plot(Phi,Ys_orignal,'--b')


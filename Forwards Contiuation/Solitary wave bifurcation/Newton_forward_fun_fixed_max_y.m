function [Phi,Xs,Ys,Xt,Yt,Theta,Theta_bottom,F,P] = Newton_forward_fun_fixed_max_y(L,N,max_y,Theta_init_guess,Theta_bottom_init_guess,F,Yt_true,P,ijac)


%Set up initial guess shapes
Theta=reshape(Theta_init_guess,N,1);
Theta_bottom=reshape(Theta_bottom_init_guess,N,1);
Yt_true=reshape(Yt_true,N,1);



%%%%%%%% initialisation
Theta_mid = zeros(N-1,1);
Tau_mid = zeros(N-1,1);
Y = zeros(N,1); 
Y_mid = zeros(N-1,1);
%%%%%%%%%%


Phi = linspace(-L,L,N);
del_Phi = 2*L/(N-1);

for i=1:N-1
    Phi_mid(i) = 0.5*(Phi(i)+Phi(i+1));
end

e_pi_phi = exp(pi*Phi);  %%% Calculate this now to save computation later



for i=1:N-1
    P_mid(i) = 0.5*(P(i)+P(i+1));
end




tol = 0.000001;
eps = 0.0001;  % For numerical differentiation
err = 2*tol;



% Create vector of unknowns



w = [Theta(2:N);Theta_bottom(2:N);F];

M = length(w);

icount = 1;
while err > tol

    if icount==1
        
        jac = zeros(M,M);

        tic
        H=eqns(w);
        for i_jac=1:M
            w_pert=w;
            w_pert(i_jac)=w_pert(i_jac)+eps;
            H_pert=eqns(w_pert);

            jac(:,i_jac)=(H_pert-H)/eps;
        end
        toc

        jac;

        % Solve linear system J.h = -H


        %h = linsolve(jac,-transpose(H));
        h = linsolve(jac,-(H));

        
    else
        H = eqns(w);
        tic
        %h = linsolve(jac,-transpose(H));
        h = linsolve(jac,-(H));
        toc
        disp('Jacobian not computed!')
    end
    
    w = w + h;
    
    err = norm(H)
    
    if ijac~=1
        icount = icount + 1;
    end
end

% Rerun the function
% with the converged w

H = eqns(w);

Theta = [0;w(1:N-1)];
Theta_bottom=[0;w(N:end-1)];
F=w(end);

% Surface shape

Xs(N) = L;
for i = 1:N-1
    Xs(N-i) = Xs(N-i+1) - (exp(-Tau_mid(N-i))*cos(Theta_mid(N-i))*del_Phi);
end

Ys(N) = Y(N);
for i = 1 : N - 1
    Ys(N-i) = Ys(N-i+1) - (exp (-Tau_mid(N-i)) * sin(Theta_mid(N-i)) * del_Phi );    %Tau_mid must be calculated first in previous loop.
end

% Topography shape

for i=1:N-1
    Tau_integrand_surface = 1./(1+exp(pi*Phi_mid(i))./e_pi_phi);
    Tau_integrand_bottom = 1./(-1+exp(pi*Phi_mid(i))./e_pi_phi);
    Tau_mid_t(i) = trapz(Phi,(transpose(Theta).*Tau_integrand_surface)) + trapz(Phi,(transpose(Theta_bottom).*Tau_integrand_bottom))  ;
end

for i=1:N-1
    Theta_mid_t(i) = 0.5*(Theta_bottom(i) + Theta_bottom(i+1));
end

Xt(N) = L;
for i = 1:N-1
    Xt(N-i) = Xt(N-i+1) - (exp(-Tau_mid_t(N-i))*cos(Theta_mid_t(N-i))*del_Phi);
end


Yt(N) = 0;
for i = 1 : N - 1
    Yt(N-i) = Yt(N-i+1) - (exp (-Tau_mid_t(N-i)) * sin(Theta_mid_t(N-i)) * del_Phi );    %Tau_mid must be calculated first in previous loop.
end








function G = eqns(w)


    Theta = [0;w(1:N-1)];
    Theta_bottom=[0;w(N:end-1)];
    F=w(end);


    Theta_mid = zeros(N-1,1);
    Tau_mid = zeros(N-1,1);
    Y = zeros(N,1); Y(N) = 1;
    Y_mid = zeros(N-1,1);
    Yt = zeros(N,1); Yt(N) = 0;


    for i = 1:N-1
        Theta_mid(i) = 0.5*(Theta(i)+Theta(i+1));
    end

    for i = 1 : N - 1
        Tau_integrand_surface = 1./(1-exp(pi*Phi_mid(i))./e_pi_phi);
        Tau_integrand_bottom = 1./(1+exp(pi*Phi_mid(i))./e_pi_phi);
        Tau_mid(i) = trapz(Phi,(transpose(Theta).*Tau_integrand_surface)) - trapz(Phi,(transpose(Theta_bottom).*Tau_integrand_bottom))  ;
    end

    for i = 1 : N - 1
        Y(N-i) = Y(N-i+1) - (exp (-Tau_mid(N-i)) * sin(Theta_mid(N-i)) * del_Phi );    %Tau_mid must be calculated first in previous loop.
    end

    for i = 1 : N - 1
        Y_mid(i) = 0.5*(Y(i) + Y(i+1)); %Y must first be calculated in previous loop
    end



    % Topography shape
    for i=1:N-1
        Tau_integrand_surface = 1./(1+exp(pi*Phi_mid(i))./e_pi_phi);
        Tau_integrand_bottom = 1./(-1+exp(pi*Phi_mid(i))./e_pi_phi);
        Tau_mid_t(i) = trapz(Phi,(transpose(Theta).*Tau_integrand_surface)) + trapz(Phi,(transpose(Theta_bottom).*Tau_integrand_bottom))  ;
    end

    for i=1:N-1
        Theta_mid_t(i) = 0.5*(Theta_bottom(i) + Theta_bottom(i+1));
    end

    for i = 1 : N - 1
        Yt(N-i) = Yt(N-i+1) - (exp (-Tau_mid_t(N-i)) * sin(Theta_mid_t(N-i)) * del_Phi );    %Tau_mid must be calculated first in previous loop.
    end





    %%% equations


    G=zeros(2*N-1,1);
    for j = 1 : N - 1
        G(j) = exp(2*Tau_mid(j)) + 2/(F^2) * (Y_mid(j) + P_mid(j)) - 1 - 2/F^2;
    end


    G(N:end-1)=Yt(1:end-1)- Yt_true(1:end-1);
    Nz = (N+1)/2;
    G(end) = Y(Nz) - max_y;
    
end




end

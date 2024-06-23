function [Yt,Ys,Xt,Xs] = variables_after_SVD(L,N_t,N_s,Theta_b,Theta)
%YT_AFTER_SVD_FUN Summary of this function goes here
%   Detailed explanation goes here

    %Make Thetas row vectors so trapz doesn't cry about column*row
    Theta_b=reshape(Theta_b,1,N_t);
    Theta=reshape(Theta,1,N_s);


    % Define meshes
    Phi_s = linspace(-L,L,N_s);
    del_Phi_s = 2*L/(N_s-1);
    Phi_t = linspace(-L,L,N_t);
    del_Phi_t = 2*L/(N_t-1);
    
    
    %Create midpoint vectors
    Phi_s_mid = 0.5*(Phi_s(1:N_s-1)+Phi_s(2:N_s));
    Phi_t_mid = 0.5*(Phi_t(1:N_t-1)+Phi_t(2:N_t));
    Theta_mid = 0.5*(Theta(1:N_s-1)+Theta(2:N_s));

    % Calculate this now to save computation later
    e_pi_phi_s = exp(pi*Phi_s);  
    e_pi_phi_t = exp(pi*Phi_t);




%% Bottom


    Theta_mid_t=0.5*(Theta_b(1:end-1)+Theta_b(2:end));


    Yt=zeros(N_t,1); % Initialise Yt
    Tau_mid_t=zeros(N_t,1); %Initialise Tau_mid
    Xt=zeros(N_t,1);
    
    
    %Evaluate Tau_mid
    for i=1:N_t-1
        Tau_integrand_surface = 1./(1+exp(pi*Phi_t_mid(i))./e_pi_phi_s);
        Tau_integrand_bottom = 1./(-1+exp(pi*Phi_t_mid(i))./e_pi_phi_t);
        Tau_mid_t(i) = trapz(Phi_s,(Theta.*Tau_integrand_surface)) + trapz(Phi_t,(Theta_b.*Tau_integrand_bottom))  ;
    end
    
    Yt(N_t) = 0; %BC on Yt
    Xt(N_t) = L;
    
   
  

    
    
    for i = 1 : N_t - 1
        Yt(N_t-i) = Yt(N_t-i+1) - (exp (-Tau_mid_t(N_t-i)) * sin(Theta_mid_t(N_t-i)) * del_Phi_t );    %Tau_mid must be calculated first in previous loop.
        Xt(N_t-i) = Xt(N_t-i+1) - (exp(-Tau_mid_t(N_t-i))*cos(Theta_mid_t(N_t-i))*del_Phi_t);
    end    

    % Turn Yt into a nice column vector
    Yt=reshape(Yt,N_t,1);

    
    
  %% Top
  
    Theta_mid=0.5*(Theta(1:end-1)+Theta(2:end));


    Ys=zeros(N_s,1); % Initialise Ys
    Ys(N_s) = 1; %BC on Ys
    
    Xs=zeros(N_s,1);
    
    Tau_mid=zeros(N_s,1); %Initialise Tau_mid
  
  
  
    %Evaluate Tau_mid
    for i = 1 : N_s - 1
        Tau_integrand_surface = 1./(1-exp(pi*Phi_s_mid(i))./e_pi_phi_s);
        Tau_integrand_bottom  = 1./(1+exp(pi*Phi_s_mid(i))./e_pi_phi_t);
        Tau_mid(i)            = trapz(Phi_s,(Theta.*Tau_integrand_surface))...
            - trapz(Phi_t,(Theta_b.*Tau_integrand_bottom));
    end
    
    
    
  
  
    Xs(N_s) = L;
  
    for i = 1 : N_s - 1
        Ys(N_s-i) = Ys(N_s-i+1) - (exp (-Tau_mid(N_s-i)) * sin(Theta_mid(N_s-i)) * del_Phi_s );    %Tau_mid must be calculated first in previous loop.
        Xs(N_s-i) = Xs(N_s-i+1) - (exp(-Tau_mid(N_s-i))*cos(Theta_mid(N_s-i))*del_Phi_s);
    end
  
    % Turn Ys into a nice column vector
    Ys=reshape(Ys,N_s,1);
    

    

  
  
  
  
    
    
    
    
    
    
end


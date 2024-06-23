function [Phi_s,Phi_t,Xs,Ys,Xt,Yt,Theta,P] = FUN_Inverse_Pressure_NO_TOPOGRAPHY(L,N_s,N_t,F,Theta)

    P=zeros(N_s,1);
    Theta_bottom=zeros(1,N_t);
   
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
    
    
    
    
    % Calculate tau_f midpoints
    for i = 1 : N_s - 1
        Tau_integrand_surface = 1./(1-exp(pi*Phi_s_mid(i))./e_pi_phi_s);
        Tau_integrand_bottom  = 1./(1+exp(pi*Phi_s_mid(i))./e_pi_phi_t);
        
        Tau_mid(i)            = trapz(Phi_s,(Theta'.*Tau_integrand_surface))...
            - trapz(Phi_t,(Theta_bottom.*Tau_integrand_bottom));
    end
    
    
    % Calculate y and y midpoints
    Y(N_s)      = 1;
    
    for i = 1 : N_s - 1
        Y(N_s-i) = Y(N_s-i+1) - (exp (-Tau_mid(N_s-i)) * sin(Theta_mid(N_s-i)) * del_Phi_s );    %Tau_mid must be calculated first in previous loop.
    end
    
    for i = 1 : N_s - 1
        Y_mid(i) = 0.5*(Y(i) + Y(i+1)); %Y must first be calculated in previous loop
    end
    
    
    %Calculate P at midpoints from Bernouilli eq
    P_mid=  1-Y_mid + ((F^2)/2)*(1-exp(2*Tau_mid));

    %Reconstruct Pressure
    P(1)=0;
    for i=1:N_s-1
        P(i+1)=2*P_mid(i) - P(i);
    end
    
    
    
    
    % After convergence find physical variables again
    [Xs,Xt,Ys,Yt] = variable_evaluate();
    

    
    function [Xs,Xt,Ys,Yt] = variable_evaluate()
        
        %% Surface
        Xs = zeros(N_s,1);
        Ys = zeros(N_s,1);
        Tau_mid   = zeros(N_s-1,1);

        Ys(N_s)      = 1;

        for i = 1 : N_s - 1
             Tau_integrand_surface = 1./(1-exp(pi*Phi_s_mid(i))./e_pi_phi_s);
             Tau_integrand_bottom  = 1./(1+exp(pi*Phi_s_mid(i))./e_pi_phi_t);
             Tau_mid(i)            = trapz(Phi_s,(Theta'.*Tau_integrand_surface))...
                                   - trapz(Phi_t,(Theta_bottom.*Tau_integrand_bottom));
        end

        for i = 1 : N_s - 1
            Ys(N_s-i) = Ys(N_s-i+1) - (exp (-Tau_mid(N_s-i)) * sin(Theta_mid(N_s-i)) * del_Phi_s );    %Tau_mid must be calculated first in previous loop.
            Xs(N_s-i) = Xs(N_s-i+1) - (exp(-Tau_mid(N_s-i))*cos(Theta_mid(N_s-i))*del_Phi_s); 
        end
        
        Xs = Xs - 0.5*Xs(1);

        %% Topography
        
        Xt = zeros(N_t,1);
        Yt = zeros(N_t,1);
        Tau_mid_t   = zeros(N_t-1,1);
        
        
        for i=1:N_t-1
            Tau_integrand_surface = 1./(1+exp(pi*Phi_t_mid(i))./e_pi_phi_s);
            Tau_integrand_bottom = 1./(-1+exp(pi*Phi_t_mid(i))./e_pi_phi_t);
            Tau_mid_t(i) = trapz(Phi_s,(Theta'.*Tau_integrand_surface)) + trapz(Phi_t,(Theta_bottom.*Tau_integrand_bottom))  ;
        end
    

        Theta_mid_t = 0.5*(Theta_bottom(1:end-1) + Theta_bottom(2:end));

    
        
        for i = 1:N_t-1
            Xt(N_t-i) = Xt(N_t-i+1) - (exp(-Tau_mid_t(N_t-i))*cos(Theta_mid_t(N_t-i))*del_Phi_t); 
            Yt(N_t-i) = Yt(N_t-i+1) - (exp (-Tau_mid_t(N_t-i)) * sin(Theta_mid_t(N_t-i)) * del_Phi_t );    %Tau_mid must be calculated first in previous loop.
        end
        Xt = Xt - 0.5*Xt(1);
    

       
  
    
    
    end




    function G = eqns(Theta_bottom)  %update y values and form system of equations for Newton solver

%           Old way of doing BCS, gave zero row in JAC
%         if F<1
%             Theta_bottom(end)=0;  % The 0 is the value of theta at last grid point (BC for subcritical)
%         else
%             Theta_bottom(1)=0;  % The 0 is the value of theta at last grid point (BC for supercritical)
%         end
% 
%         Tau_mid   = zeros(N_s-1,1);
%         Y         = zeros(N_s,1); 
%         Y(N_s)      = 1;
%         Y_mid     = zeros(N_s-1,1);


% 
%         
%         for j = 1 : N_s - 1
%             G(j) = exp(2*Tau_mid(j)) + 2/(F^2) * (Y_mid(j) + P_mid(j)) - 1 - 2/F^2;
%         end
%         
%         if F<1
%             G(N_s)=Theta_bottom(end);  % The 0 is the value of theta at last grid point (BC for subcritical)
%         else
%             G(N_s)=Theta_bottom(1);  % The 0 is the value of theta at last grid point (BC for supercritical)
%         end

 
    end
end
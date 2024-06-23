function [Phi,Xs,Ys,Xt,Yt,Theta,Theta_bottom,P] = surfacefinder_topography_fixedFP(L,N,F,ijac,P,Theta_bottom,Theta)
    

    Theta=reshape(Theta,N,1);



    Phi = linspace(-L,L,N);
    del_Phi = 2*L/(N-1);
    
    tol = 1E-6;
    eps = 0.0001;  % For numerical differentiation
    err = 2*tol;

    
    


        
    
    for i=1:N-1
        Phi_mid(i) = 0.5*(Phi(i)+Phi(i+1));
    end
    
    e_pi_phi = exp(pi*Phi);  %%% Calculate this now to save computation later
   
    
    for i=1:N-1
        P_mid(i) = 0.5*(P(i)+P(i+1));
    end
    
    


    
    
    
    Theta_mid=[];
    Tau_mid=[];
    Y=[];
    
    
    %%%%%%%%% compute Jacobian numerically and Newton solve
    % Create vector of unknowns
    if F>1
        w = [Theta(2:N)];
    elseif F<1
        w=[Theta(1:N-1)];
    end
    
    M = length(w);
    
    icount = 1;


    
    while err > tol

        if icount==1  
            jac = zeros(M,M);
            tic %%%how long does the Jacobian take to fill
            H=eqns(w);
            for i_jac=1:M
                w_pert=w;
                w_pert(i_jac)=w_pert(i_jac)+eps;
                H_pert=eqns(w_pert);
                
                jac(:,i_jac)=(H_pert-H)/eps;
            end
            toc
            
            jac;
                        
            %%%% Solve linear system J.h = -H
            
            tic %%%how long does solving matrix equation take

            h = linsolve(jac,-transpose(H));
            toc

        else
            H = eqns(w);
            tic
            h = linsolve(jac,-transpose(H));
            toc
            disp('Jacobian not computed!')
        end

        w = w + h;

        err = norm(H)
        
        if ijac~=1
            icount = icount + 1;    %%%this line means the Jacobian is only calculated once, on the first loop
        end
    end
    
    %%% Rerun the function with the converged w
    H = eqns(w);
    if F<1
        Theta = [w; 0];  % The 0 is the value of theta at last grid point (BC for subcritical)
    elseif F>1
        Theta = [0; w];  % The 0 is the value of theta at last grid point (BC for supercritical)
    end

    % Surface shape

    Xs(N) = 0;
    for i = 1:N-1
        Xs(N-i) = Xs(N-i+1) - (exp(-Tau_mid(N-i))*cos(Theta_mid(N-i))*del_Phi); 
    end
    Xs = Xs - 0.5*Xs(1);
    
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
    
    Xt(N) = 0;
    for i = 1:N-1
        Xt(N-i) = Xt(N-i+1) - (exp(-Tau_mid_t(N-i))*cos(Theta_mid_t(N-i))*del_Phi); 
    end
    Xt = Xt - 0.5*Xt(1);
    
    Yt(N) = 0;
    for i = 1 : N - 1
            Yt(N-i) = Yt(N-i+1) - (exp (-Tau_mid_t(N-i)) * sin(Theta_mid_t(N-i)) * del_Phi );    %Tau_mid must be calculated first in previous loop.
    end

    



    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function G = eqns(w)  %update y values and form system of equations for Newton solver


        if F<1
            Theta = [w; 0];  % The 0 is the value of theta at last grid point (BC for subcritical)
        elseif F>1
            Theta = [0; w];  % The 0 is the value of theta at last grid point (BC for supercritical)
        end


        Theta_mid = zeros(N-1,1);
        Tau_mid = zeros(N-1,1);
        Y = zeros(N,1); Y(N) = 1;
        Y_mid = zeros(N-1,1);

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
        
        for j = 1 : N - 1
            G(j) = exp(2*Tau_mid(j)) + 2/(F^2) * (Y_mid(j) + P_mid(j)) - 1 - 2/F^2;
        end
 
    end
    
end
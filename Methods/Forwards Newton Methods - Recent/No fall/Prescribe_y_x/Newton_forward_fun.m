function [Phi,Xs,Ys,Xt,Yt,Theta,Theta_bottom,F] = Newton_forward_fun(L,N,F,amplitude,b_width,ijac,P,Theta_init_guess,Theta_b_init_guess)


    %Set up initial guess for theta
    Theta=reshape(Theta_init_guess,N,1);
    Theta_bottom=reshape(Theta_b_init_guess,N,1);
    
    Phi = linspace(-L,L,N);
    del_Phi = 2*L/(N-1);
    
    tol = 1E-6;
    eps = 0.0001;  % For numerical differentiation
    err = 2*tol;

    Theta_mid=[];
    Tau_mid=[];
    Y=[];
    
       
    
    Phi_mid=0.5*(Phi(1:end-1)+Phi(2:end));
    
    
    e_pi_phi = exp(pi*Phi);  %%% Calculate this now to save computation later
   
    
    P_mid=0.5*(P(1:end-1)+P(2:end));

    
    




    
    %%%%%%%%% compute Jacobian numerically and Newton solve
    % Create vector of unknowns
    if F>1
        w = [Theta(2:N);Theta_bottom(2:N)];
    else
        w=[Theta(1:N-1);Theta_bottom(1:N-1)];
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
            
%             jac1=jac;
%             save('JAC1.mat','jac1')
            
            %jac;
                        
            %%%% Solve linear system J.h = -H
            
            tic %%%how long does solving matrix equation take

            
            %h = linsolve(jac,-transpose(H));
            h = linsolve(jac,-(H));
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
        Theta = [w(1:N-1); 0];  % The 0 is the value of theta at last grid point (BC for subcritical)
        Theta_bottom=[w(N:end); 0];
    elseif F>1
        Theta = [0; w(1:N-1)];  % The 0 is the value of theta at last grid point (BC for supercritical)
        Theta_bottom=[0; w(N:end)];
    end
    
    
    

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

    
    
   


    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function G = eqns(w)  %update y values and form system of equations for Newton solver

        
        if F<1
            Theta = [w(1:N-1); 0];  % The 0 is the value of theta at last grid point (BC for subcritical)
            Theta_bottom=[w(N:end);0];
        elseif F>1
            Theta = [0; w(1:N-1)];  % The 0 is the value of theta at last grid point (BC for supercritical)
            Theta_bottom=[0;w(N:end)];
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

        
    dYdx_Theta_b=atan(-2*amplitude*(b_width^2)*Xt.*exp(-(b_width*Xt).^2))';
        
        
        
        
        
        G=zeros(2*N-2,1);
        
        
        for j = 1 : N - 1
            G(j) = exp(2*Tau_mid(j)) + 2/(F^2) * (Y_mid(j) + P_mid(j)) - 1 - 2/F^2;
        end
        
        
        %G(N:end)=Theta_bottom(2:end)-dYdx_Theta_b(2:end);
        
         if F<1
            G(N:end)=Theta_bottom(1:end-1)-dYdx_Theta_b(1:end-1);
         else
             G(N:end)=Theta_bottom(2:end)-dYdx_Theta_b(2:end);
         end

 
    end
    
end
function [Phi,Xs,Ys,Xt,Yt,Theta,Theta_bottom,P] = surfacefinder_fixedFP(L,N,F,iload,isave,ijac,Pressure,Topography)
    
    Phi = linspace(-L,L,N);
    del_Phi = 2*L/(N-1);
    
    tol = 1E-6;
    eps = 0.0001;  % For numerical differentiation
    err = 2*tol;

    
    
    if Topography ==1
        Theta_bottom=-0.01*Phi.*exp(-(Phi).^2);
    else
        Theta_bottom=zeros(1,N);
    end
    
    if Pressure ==1
        P = 0.01*exp(-(Phi).^2); % initial pressure distribution
    else 
        P = zeros(N,1);
    end
    
        
    
    for i=1:N-1
        Phi_mid(i) = 0.5*(Phi(i)+Phi(i+1));
    end
    
    e_pi_phi = exp(pi*Phi);  %%% Calculate this now to save computation later
   
    
    for i=1:N-1
        P_mid(i) = 0.5*(P(i)+P(i+1));
    end
    
    
    if iload==0
        Theta_init = zeros(N,1); %Initial guess for theta
        Theta = Theta_init; %Initialises Theta
    end
    if iload==1
        Ntemp   = N;
        Phitemp = Phi;
        Theta   = [];
        Phi     = [];
        load saved.mat
        if Ntemp>N
            Thetatemp = spline(Phi,Theta,Phitemp);
            %plot(Phi,Theta,'or',Phitemp,Theta_new,'-xb')
            Theta = transpose(Thetatemp);
            Phi   = transpose(Phitemp);
        end    
        N = Ntemp;
    end
    if iload==2
        load Forunforced_thetaguess.mat
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
            for i_jac=1:M
                H = eqns(w);
                for j_jac=1:M
                    w(j_jac) = w(j_jac) + eps;
                    Hp   = eqns(w);
                    jac(i_jac,j_jac) = (Hp(i_jac)-H(i_jac))/eps;
                    w(j_jac) = w(j_jac) - eps;
                end
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
        Tau_mid_t(i) = trapz(Phi,(transpose(Theta).*Tau_integrand_surface)) + trapz(Phi,(Theta_bottom.*Tau_integrand_bottom))  ;
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

    
    
    %%%save the data
    if isave==1 %%%saves data to be fed back into code from iload
     save('saved2.mat','Theta','N','Phi','P','F','X','Y')
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
            Tau_mid(i) = trapz(Phi,(transpose(Theta).*Tau_integrand_surface)) - trapz(Phi,(Theta_bottom.*Tau_integrand_bottom))  ;

    %     % Check trapezium rule
    %     
    %         sum = 0.5*del_Phi*Theta(1)*Tau_integrand(1);
    %         for iq=2:N-1
    %             sum = sum + del_Phi*Theta(iq)*Tau_integrand(iq);  
    %         end
    %         sum = sum + 0.5*del_Phi*Theta(N)*Tau_integrand(N);
    % 
    %         if abs(sum)>0.001
    %             format long
    %             qcheck1 = sum
    %             qcheck2 = Tau_mid(i)
    %             qcheck3 = abs(sum-Tau_mid(i))
    %             format short
    %             pause
    %         end


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
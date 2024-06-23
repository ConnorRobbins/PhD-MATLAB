function [Phi,X,Y,P,Theta] = pressurefinder_fixedFTheta(L,N,F,iload,isave,Pzero_origin)

    tol = 1E-11;
    eps = 0.0001;  % For numerical differentiation
    err = 2*tol;

    %%%%
    % Grid
    %%%%

    Phi = linspace(-L,L,N);
    del_Phi = 2*L/(N-1);

    for i=1:N-1
        Phi_mid(i) = 0.5*(Phi(i)+Phi(i+1));
    end

    e_pi_phi = exp(pi*Phi);  %%% Calculate this now to save computation later

    if iload==0
        Theta_init = zeros(N,1); %Initial guess for theta
        Theta = Theta_init; %Initialises Theta

        %%%%%%%% Try a Gaussian-ish free surface
        Theta = -2*Phi.*exp(-Phi.^2);
        Theta = transpose(Theta);
    end


    Phi=transpose(Phi);

    %%%%%%%%%%%%%%%% Initial guess Pressure distribution
    P = zeros(N,1);     
    %P = (0.01*exp(-(Phi+5).^2)); 
    %P = 0.01*exp(-Phi.^2);
    %%%%


    for i=1:N-1
        Theta_mid(i) = 0.5*(Theta(i)+Theta(i+1));
    end


    Tau_mid = zeros(N-1,1); %%%initialise Tau_mid
    for i = 1 : N - 1  %%% Fill in Tau_mid
    %       Tau_integrand = e_pi_phi./(e_pi_phi-exp(pi*Phi_mid(i)));
        Tau_integrand = 1./(1-exp(pi*Phi_mid(i))./e_pi_phi);
        Tau_mid(i) = trapz(Phi,(transpose(Theta).*Tau_integrand));  
    end
    
    
    Y = zeros(N,1); Y(N) = 1;   %%%initialise Y and set downstream boundary condition
    Y_mid = zeros(N-1,1);

 
    for i = 1 : N - 1
        Y(N-i) = Y(N-i+1) - (exp (-Tau_mid(N-i)) * sin(Theta_mid(N-i)) * del_Phi );    %Tau_mid must be calculated first in previous loop.
    end

    for i = 1 : N - 1
        Y_mid(i) = 0.5*(Y(i) + Y(i+1)); %Y must first be calculated in previous loop
    end
    



    % Create vector of unknowns
    w = [P(1:N-1)];
    M = length(w);


    Pinit=P;   %%%store inital guess
    
    icount = 1;

    
    %%%%%%% Newton Solver 

    while err > tol

        if icount==1    %evaluate and fill in Jacobian

            jac = zeros(M,M);
            tic
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

            % Solve linear system J.h = -H
            tic
            h = linsolve(jac,-transpose(H));
            toc

        else   %%%% Use previously calculated Jacobian
            H = eqns(w);
            tic
            h = linsolve(jac,-transpose(H));
            toc
            disp('Jacobian not computed!')
        end

        w = w + h;    
        err = norm(H)    
        icount = icount + 1;
    end
    

    % Rerun the function with the converged w
    H = eqns(w);
    P = [w; 0];

    X(N) = 0;
    for i = 1:N-1
        X(N-i) = X(N-i+1) - (exp(-Tau_mid(N-i))*cos(Theta_mid(N-i))*del_Phi); 
    end
    X = X - 0.5*X(1);



    if isave==1 %%%saves data to be fed back into code from iload
        save('saved2.mat','Theta','N','Phi')
    end


    function G = eqns(w)

        P = [w; 0];  % The 0 is the value of theta at last grid point (BC for subcritical)
        
        
        if Pzero_origin==1
            if mod(N,2)==0
                P(N/2)=0;
                P((N/2)+1)=0;
            else
                P((N+1)/2)=0;
            end
        end
        
        
        for i=1:N-1
            P_mid(i) = 0.5*(P(i)+P(i+1));
        end


        for j = 1 : N - 1
            G(j) = exp(2*Tau_mid(j)) + 2/(F^2) * (Y_mid(j) + P_mid(j)) - 1 - 2/F^2;
        end

    end
end
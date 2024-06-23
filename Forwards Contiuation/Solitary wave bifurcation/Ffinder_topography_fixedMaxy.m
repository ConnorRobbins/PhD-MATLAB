function [Phi,X,Y,F,Theta,P] = Ffinder_fixedMaxyP(L,N,iload,isave,mxy,iF,Theta,F)
 iF




% Newton Solver    
%



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

% Pressure distribution

%P = 0.01*exp(-Phi.^2); zeros(N,1); % initial pressure distribution
P = zeros(N,1);

for i=1:N-1
    P_mid(i) = 0.5*(P(i)+P(i+1));
end


Theta_bottom=-2*0.02*Phi.*exp(-(Phi).^2);

tol = 0.000001;
eps = 0.0001;  % For numerical differentiation
err = 2*tol;



% Create vector of unknowns

w = [Theta(1:N-1); F];
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

        tic
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
    
%    icount = icount + 1;
end

% Rerun the function
% with the converged w

H = eqns(w);

Theta = [w(1:end-1); 0];

X(N) = 0;
for i = 1:N-1
    X(N-i) = X(N-i+1) - (exp(-Tau_mid(N-i))*cos(Theta_mid(N-i))*del_Phi); 
end
X = X - 0.5*X(1);

% subplot(2,1,1)
% plot(Phi,Y,'-ob')
% ylabel('y')
% subplot(2,1,2)
% plot(Phi,Theta,'-xr')
% xlabel('Phi')
% ylabel('theta')
%  





% save('bifdig','Fr','mxx')
% 
% load bifdig-1
% plot(Fr,mxx,'-xr')

% string1   = 'output';
% FilenameN = strcat(string1,num2str(N));
% save('saved.mat','Phi','X','Y','N')
% save(FilenameN,'Phi','X','Y','Theta','N') %%%for saving outputs for later
% analysis

        if F<1
            Theta = [w; 0];  % The 0 is the value of theta at last grid point (BC for subcritical)
        elseif F>1
            Theta = [0; w];  % The 0 is the value of theta at last grid point (BC for supercritical)
        end
%
%
%




if isave==1 %%%saves data to be fed back into code from iload
    save('saved2.mat','Theta','N','Phi')
end


function G = eqns(w)

    
    F = w(end);
    if F<1
        Theta = [w(1:end-1); 0];  % The 0 is the value of theta at last grid point (BC for subcritical)
    elseif F>1
        Theta = [0; w(1:end-1)];  % The 0 is the value of theta at last grid point (BC for supercritical)
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
             Tau_mid(i) = trapz(Phi,(transpose(Theta).*Tau_integrand_surface)) - trapz(Phi,(Theta_bottom.*Tau_integrand_bottom));
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

    Nz = (N+1)/2;
    G(N) = Y(Nz) - mxy;
    
end
end

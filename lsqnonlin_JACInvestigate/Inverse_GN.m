%close all
clear
clc

L=15;
N_s_points = 41;
N_t_points = 41;

Froude=1.1;
amplitude=0.02;
b_width=sqrt(1); %sqrt(what Ben had in previous papers)


FreeSurfaceType=1;
Pressure=0;



%Initial guess
x0=zeros(1,N_t_points);

[Phi_s,Phi_t,Xs,Ys,Xt,Yt,Theta,Theta_bottom,P,jac_inv_GN]= funcaller(L,N_s_points,N_t_points,Froude,x0,FreeSurfaceType,Pressure,amplitude,b_width);

jac_inv_GN=full(jac_inv_GN);

figure
%% Plotting

titlestring=strcat(' F=',num2str(Froude),';',' Ns=',num2str(N_s_points),';',' Nt=',num2str(N_t_points),';',' a=',num2str(amplitude),';',' b=',num2str(b_width))
plot(Phi_t,Yt,'-r')
hold on
plot(Phi_s,Ys-1,'-b')
%plot(Phi_s,Ys,'-r')
legend('Topography','Surface - 1')
suptitle(strcat('Inverse GN:',titlestring))


%%


function [Phi_s,Phi_t,Xs,Ys,Xt,Yt,Theta,Theta_bottom,P,jac_inv_GN] = funcaller(L,N_s,N_t,F,x0,FreeSurfaceType,Pressure,amplitude,b_width)

    P=zeros(N_s,1);
   
% Define meshes
    Phi_s = linspace(-L,L,N_s);
    del_Phi_s = 2*L/(N_s-1);
    
    Phi_t = linspace(-L,L,N_t);
    del_Phi_t = 2*L/(N_t-1);
    
    
% Define variables on meshes
    if FreeSurfaceType==0
        Theta=zeros(1,N_s); %flat
    end
    if FreeSurfaceType==1
        Theta = atan(-2*amplitude*b_width*Phi_s.*exp(-(b_width*Phi_s).^2));  % Gaussian
    end    
    if FreeSurfaceType==2
        Theta = atan(-amplitude*b_width*2*sech(b_width*Phi_s).^2.*tanh(b_width*Phi_s)); %sech^2
    end  
    
    if Pressure ==1
        P = 0.01*exp(-(Phi_s).^2); % initial pressure distribution
    else 
        P = zeros(N_s,1);
    end
    
    
    
    %Create midpoint vectors
    Phi_s_mid = 0.5*(Phi_s(1:N_s-1)+Phi_s(2:N_s));
    Phi_t_mid = 0.5*(Phi_t(1:N_t-1)+Phi_t(2:N_t));
    P_mid = 0.5*(P(1:N_s-1)+P(2:N_s));
    Theta_mid = 0.5*(Theta(1:N_s-1)+Theta(2:N_s));
    
    % Calculate this now to save computation later
    e_pi_phi_s = exp(pi*Phi_s);  
    e_pi_phi_t = exp(pi*Phi_t);
    
    
    
    [topo_vec,resnorm,res,eflag,output1,lambda,jac_inv_GN]=lsqnonlin(@eqns,x0);
    
    resnorm
    res;
    eflag
    output1
    
    Theta_bottom=topo_vec;
    
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
             Tau_mid(i)            = trapz(Phi_s,(Theta.*Tau_integrand_surface))...
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
            Tau_mid_t(i) = trapz(Phi_s,(Theta.*Tau_integrand_surface)) + trapz(Phi_t,(Theta_bottom.*Tau_integrand_bottom))  ;
        end
    

        Theta_mid_t = 0.5*(Theta_bottom(1:end-1) + Theta_bottom(2:end));

    
        
        for i = 1:N_t-1
            Xt(N_t-i) = Xt(N_t-i+1) - (exp(-Tau_mid_t(N_t-i))*cos(Theta_mid_t(N_t-i))*del_Phi_t); 
            Yt(N_t-i) = Yt(N_t-i+1) - (exp (-Tau_mid_t(N_t-i)) * sin(Theta_mid_t(N_t-i)) * del_Phi_t );    %Tau_mid must be calculated first in previous loop.
        end
        Xt = Xt - 0.5*Xt(1);
    

       
  
    
    
    end




    function G = eqns(Theta_bottom)  %update y values and form system of equations for Newton solver
% 
%         if F<1
%             Theta_bottom(end)=0;  % The 0 is the value of theta at last grid point (BC for subcritical)
%         else
%             Theta_bottom(1)=0;  % The 0 is the value of theta at last grid point (BC for supercritical)
%         end

        Tau_mid   = zeros(N_s-1,1);
        Y         = zeros(N_s,1); 
        Y(N_s)      = 1;
        Y_mid     = zeros(N_s-1,1);

        for i = 1 : N_s - 1
             Tau_integrand_surface = 1./(1-exp(pi*Phi_s_mid(i))./e_pi_phi_s);
             Tau_integrand_bottom  = 1./(1+exp(pi*Phi_s_mid(i))./e_pi_phi_t);
             
             
             
             Tau_mid(i)            = trapz(Phi_s,(Theta.*Tau_integrand_surface))...
                                   - trapz(Phi_t,(Theta_bottom.*Tau_integrand_bottom));
        end

        for i = 1 : N_s - 1
            Y(N_s-i) = Y(N_s-i+1) - (exp (-Tau_mid(N_s-i)) * sin(Theta_mid(N_s-i)) * del_Phi_s );    %Tau_mid must be calculated first in previous loop.
        end

        for i = 1 : N_s - 1 
            Y_mid(i) = 0.5*(Y(i) + Y(i+1)); %Y must first be calculated in previous loop
        end
        
        for j = 1 : N_s - 1
            G(j) = exp(2*Tau_mid(j)) + 2/(F^2) * (Y_mid(j) + P_mid(j)) - 1 - 2/F^2;
        end
        
        if F<1
            G(N_s)=Theta_bottom(end);  % The 0 is the value of theta at last grid point (BC for subcritical)
        else
            G(N_s)=Theta_bottom(1);  % The 0 is the value of theta at last grid point (BC for supercritical)
        end
 
    end
end
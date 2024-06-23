%clear
L_width=15;
N_s_points = 41;
N_t_points = 41;

Froude=0.8;
amplitude=0.02;
b_width=sqrt(1); %sqrt(what Ben had in previous papers)


FreeSurfaceType=1;
Pressure=0;
isave=0;
Topography=1;

makeplots='no';






[Phi_s,Phi_t,Xs,Ys,Xt,Yt,Theta_surface,Theta_bottom,P,jac_GN]=forward_fun(L_width,N_s_points,N_t_points,Froude,isave,Pressure,Topography,amplitude,b_width);
jajac_GN=full(jac_GN);
%titlestring=strcat(' F=',num2str(Froude),';',' N=',num2str(N_points),';',' a=',num2str(amplitude),';',' b=',num2str(b_width));



%%   Plotting
if strcmp(makeplots,'yes')==1;

    figure 
    suptitle(strcat('Forward GN:',titlestring))
    
    if Pressure == 0;

        subplot(2,1,1)
        %plot(Xs,Ys-1,'-ob')
        plot(Phi_s,Ys-1,'-ob')
        %ylim([-0.1,0])
        hold on
        %plot(Xt,Yt,'-or')
        plot(Phi_t,Yt,'-or')
        legend('Surface-1','Topography')
        %xlabel('X')
        ylabel('$y$','interpreter','latex','FontSize',16)
        %xlim([Xt(1) Xt(end)])

        subplot(2,1,2)
        %plot(Xs,Theta_surface,'-xb')
        plot(Phi_s,Theta_surface,'-xb')
        hold on
        %plot(Xt,Theta_bottom,'-xr')
        plot(Phi_t,Theta_bottom,'-xr')
        legend('Surface','Topography')
        %title('Forward Problem Solver')
        %xlabel('X')
        xlabel('$\phi$','interpreter','latex','FontSize',16)
        ylabel('$\theta$','interpreter','latex','FontSize',16)
        %xlim([Xt(1) Xt(end)])

    end

    if Pressure == 1;

        subplot(3,1,1)
        %plot(Xs,Ys-1,'-ob')
        plot(Phi_s,Ys-1,'-ob')
        hold on
        %plot(Xt,Yt,'-or')
        plot(Phi_t,Yt,'-or')
        legend('Surface-1','Topography')
        %xlabel('X')
        ylabel('$y$','interpreter','latex','FontSize',16)
        %xlim([Xt(1) Xt(end)])

        subplot(3,1,2)
        %plot(Xs,Theta_surface,'-xb')
        plot(Phi_s,Theta_surface,'-xb')
        hold on
        %plot(Xt,Theta_bottom,'-xr')
        plot(Phi_t,Theta_bottom,'-xr')
        legend('Surface','Topography')
        %title('Forward Problem Solver')
        %xlabel('X')
        ylabel('$\theta$','interpreter','latex','FontSize',16)
        %xlim([Xt(1) Xt(end)])

        subplot(3,1,3)
        plot(Phi_s,P)
        ylabel('$P$','interpreter','latex','FontSize',16)
        xlabel('$\phi$','interpreter','latex','FontSize',16)

    end

end

%%




function [Phi_s,Phi_t,Xs,Ys,Xt,Yt,Theta,Theta_bottom,P,jac_GN] = forward_fun(L,N_s,N_t,F,isave,Pressure,Topography,amplitude,b_width)
    
    Phi_s = linspace(-L,L,N_s);
    del_Phi_s = 2*L/(N_s-1);
    
    Phi_t = linspace(-L,L,N_t);
    del_Phi_t = 2*L/(N_t-1);
    
    tol = 1E-6;
    eps = 0.0001;  % For numerical differentiation
    err = 2*tol;

    
    
    if Topography ==1
        Theta_bottom=atan(-2*amplitude*b_width*Phi_t.*exp(-(b_width*Phi_t).^2));
    else
        Theta_bottom=zeros(1,N_t);
    end
    
    if Pressure ==1
        P = -0.1*exp(-(Phi_s).^2); % initial pressure distribution
    else 
        P = zeros(N_s,1);
    end
    
        
    
    for i=1:N_s-1
        Phi_mid_s(i) = 0.5*(Phi_s(i)+Phi_s(i+1));
    end
    
    for i=1:N_t-1
        Phi_mid_t(i) = 0.5*(Phi_t(i)+Phi_t(i+1));
    end
    
    e_pi_phi_s = exp(pi*Phi_s);  %%% Calculate this now to save computation later
    e_pi_phi_t = exp(pi*Phi_t);
    
    for i=1:N_s-1
        P_mid(i) = 0.5*(P(i)+P(i+1));
    end
    
    Theta_init=zeros(N_s-1,1);
    
    %%% Call lsqnonlin
    
    [Theta_solved,resnorm,res,eflag,output1,lambda,jac_GN]=lsqnonlin(@eqns,Theta_init);
    
    Theta=Theta_solved;
    
    if F<1
        Theta = [Theta; 0];  % The 0 is the value of theta at last grid point (BC for subcritical)
    elseif F>1
        Theta = [0; Theta];  % The 0 is the value of theta at last grid point (BC for supercritical)
    end
    
    [Xs,Xt,Ys,Yt]=variable_evaluate();




    function [Xs,Xt,Ys,Yt] = variable_evaluate()
        
        %% Surface
        Xs = zeros(N_s,1);
        Ys = zeros(N_s,1);
        Tau_mid   = zeros(N_s-1,1);

        Ys(N_s)      = 1;
        
        Theta_mid = 0.5*(Theta(1:end-1) + Theta(2:end));

        for i = 1 : N_s - 1
             Tau_integrand_surface = 1./(1-exp(pi*Phi_mid_s(i))./e_pi_phi_s);
             Tau_integrand_bottom  = 1./(1+exp(pi*Phi_mid_s(i))./e_pi_phi_t);
             Tau_mid(i)            = trapz(Phi_s,((Theta.').*Tau_integrand_surface))...
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
            Tau_integrand_surface = 1./(1+exp(pi*Phi_mid_t(i))./e_pi_phi_s);
            Tau_integrand_bottom = 1./(-1+exp(pi*Phi_mid_t(i))./e_pi_phi_t);
            Tau_mid_t(i) = trapz(Phi_s,((Theta.').*Tau_integrand_surface)) + trapz(Phi_t,(Theta_bottom.*Tau_integrand_bottom))  ;
        end
    

        Theta_mid_t = 0.5*(Theta_bottom(1:end-1) + Theta_bottom(2:end));

    
        
        for i = 1:N_t-1
            Xt(N_t-i) = Xt(N_t-i+1) - (exp(-Tau_mid_t(N_t-i))*cos(Theta_mid_t(N_t-i))*del_Phi_t); 
            Yt(N_t-i) = Yt(N_t-i+1) - (exp (-Tau_mid_t(N_t-i)) * sin(Theta_mid_t(N_t-i)) * del_Phi_t );    %Tau_mid must be calculated first in previous loop.
        end
        Xt = Xt - 0.5*Xt(1);
    
    
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function G = eqns(Theta)  %update y values and form system of equations for Newton solver

        Theta_mid = zeros(N_s-1,1);
        Tau_mid = zeros(N_s-1,1);
        Y = zeros(N_s,1); 
        Y_mid = zeros(N_s-1,1);
        %Boundary conditions
        Y(N_s) = 1;
        %Theta(N_s)=0;
        
        if F<1
            Theta = [Theta; 0];  % The 0 is the value of theta at last grid point (BC for subcritical)
        elseif F>1
            Theta = [0; Theta];  % The 0 is the value of theta at last grid point (BC for supercritical)
        end
        
        
        
        

        for i = 1:N_s-1
            Theta_mid(i) = 0.5*(Theta(i)+Theta(i+1));
        end

        for i = 1 : N_s - 1
            Tau_integrand_surface = 1./(1-exp(pi*Phi_mid_s(i))./e_pi_phi_s);
            Tau_integrand_bottom = 1./(1+exp(pi*Phi_mid_s(i))./e_pi_phi_t);
            Tau_mid(i) = trapz(Phi_s,((Theta.').*Tau_integrand_surface)) - trapz(Phi_t,(Theta_bottom.*Tau_integrand_bottom))  ;
   

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
 
%         if F<1
%             G(N_s) = Theta(end);  % The 0 is the value of theta at last grid point (BC for subcritical)
%         else
%             G(N_s) = Theta(1);  % The 0 is the value of theta at last grid point (BC for supercritical)
%         end
    end
    
end

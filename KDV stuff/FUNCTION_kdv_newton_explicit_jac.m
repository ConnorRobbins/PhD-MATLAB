function [eta_solved] = FUNCTION_kdv_newton_explicit_jac(N,forcing,mu,del_x,eta_init_guess,ijac)

    tol = 1E-6;
    eps = 0.0001;  % For numerical differentiation
    err = 2*tol;



    %%%%%%%%% compute Jacobian numerically and Newton solve
    % Create vector of unknowns
    
    w=reshape(eta_init_guess,N,1);
    forcing=reshape(forcing,N,1);

    
    M = length(w);
    
    icount = 1;

    
    while err > tol

        if icount==1  
            jac = zeros(M,M);
            tic %%%how long does the Jacobian take to fill
            eta=w;
            
            %%%%%%%% System of N equations,
            %
            %   EQNS (1, N) are other conditions i.e. eta(1)=0
            %   EQNS 2:N-1 are finite difference kdv
            %
            %
            
            
            
            %fill derivatives of EQNS (1,  N)
            jac(1,1)=1; %eta(1)=0

            if mu<0
                jac(N,2)=1; %eta(2)=0
            else
                jac(N,N)=1; %eta(N)=0
            end


            
            %fill derivatives of EQNS 2-N-1
            for i=2:N-1
                %wrt eta_j
                jac(i,i)=mu-1.5*eta(i) +1/(3*(del_x^2));
                jac(i,i-1)=-1/(6*(del_x^2));
                jac(i,i+1)=-1/(6*(del_x^2));
            end
            toc
            
      
            
                        
            %%%% Solve linear system J.h = -H
            H=eqns(w);
            tic %%%how long does solving matrix equation take
            h = linsolve(jac,-(H));
            toc

        else
            H = eqns(w);
            tic
            h = linsolve(jac,-(H));
            toc
            disp('Jacobian not computed!')
        end

        w = w + h;
       
        err = norm(H)
        
        if ijac~=1
            icount = icount + 1;    %%%this line means the Jacobian is only calculated once, on the first loop
        end
    end
    
    eta_solved=w;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function H = eqns(w)  %form system of equations for Newton solver

        eta=w;
        
        H=zeros(N,1);
        
        H(1)=eta(1);
        H(2:N-1)= mu*eta(2:N-1)  - 3*(eta(2:N-1).^2)/4   - ((eta(3:N) -2*eta(2:N-1) +eta(1:N-2)) /  (6*del_x^2)) - 0.5*forcing(2:N-1);
        
        if mu<0
            H(N)=eta(2);
        else
            H(N)=eta(end);
        end
 
    end
    
end
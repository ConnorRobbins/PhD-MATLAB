function [eta,c] = FUNCTION_kdv_subcritical_trapped_waves(N,a,b,mu,x,del_x,ijac,initial_unknowns);

    tol = 1E-6;
    err = 2*tol;



    %%%%%%%%% compute Jacobian numerically and Newton solve
    % Create vector of unknowns
    
    w=reshape(initial_unknowns,N+1,1);
  

    
    M = length(w);
    
    icount = 1;

    
    while err > tol

        if icount==1  
            jac = zeros(M,M);
            tic %%%how long does the Jacobian take to fill
            eta=w(1:end-1);
            c=w(end);
            
            %%%%%%%% System of N equations,
            %
            %   EQNS (1, N) are other conditions i.e. eta(1)=0
            %   EQNS 2:N-1 are finite difference kdv
            %
            %
            
            
            
            %fill derivatives of EQNS (1,N,N+1)
            jac(1,1)=1; %eta(1)=0
            jac(N,N-1)=1; %eta(end-1)=0
            jac(N+1,N)=1;%eta(end)=0



            
            %fill derivatives of EQNS 2-N-1
            for i=2:N-1
                %wrt eta_j
                jac(i,i)=mu-1.5*eta(i) +1/(3*(del_x^2));
                jac(i,i-1)=-1/(6*(del_x^2));
                jac(i,i+1)=-1/(6*(del_x^2));
                
                %wrt to c
                jac(i,N+1)= -2*a*(b^2)*((x(i)-c)*exp( -(b*(x(i)-c))^2 ) - (x(i)+c)*exp( -(b*(x(i)+c))^2 ));
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
    
    eta=w(1:end-1);
    c=w(end);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function H = eqns(w)  %form system of equations for Newton solver

        
        eta=w(1:end-1);
        c=w(end);
        
        forcing=a*exp(-((b*(x-c)).^2)) +a*exp(-((b*(x+c)).^2));
        forcing=reshape(forcing,N,1);
        
        H=zeros(N+1,1);
        
        H(1)=eta(1);
        H(2:N-1)= mu*eta(2:N-1)  - 3*(eta(2:N-1).^2)/4   - ((eta(3:N) -2*eta(2:N-1) +eta(1:N-2)) /  (6*del_x^2)) - 0.5*forcing(2:N-1);
        
        
        H(N)=eta(end-1);
        H(N+1)=eta(end);

 
    end
    
end
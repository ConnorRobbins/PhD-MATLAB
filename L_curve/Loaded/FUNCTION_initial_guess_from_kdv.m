function [eta_solved] = FUNCTION_initial_guess_from_kdv(N,forcing,F,del_x,ijac)

    tol = 1E-6;
    err = 2*tol;
    NewtonLoopMaxIter=1E4;


    
    %%%%%%%%% compute Jacobian numerically and Newton solve
    % Create vector of unknowns
    
    eta_init_guess=zeros(1,N);
    w=reshape(eta_init_guess,N,1);
    forcing=reshape(forcing,N,1);

    
    M = length(w);
    
    icount = 1;
    NewtonLoopIter=0;
    mu=F-1;

    
    while err > tol
        NewtonLoopIter=NewtonLoopIter+1;
        
        if icount==1  
            jac = zeros(M,M);
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
      
             
            
            
        %%%% Solve linear system J.h = -H
            H=eqns(w);
            [h,recipcond] = linsolve(jac,-(H));

        else
            H = eqns(w);
            [h,recipcond] = linsolve(jac,-(H));
            %disp('Jacobian not computed!')
        end

        w = w + h;
       
        err = norm(H);

        
        if ijac~=1
            icount = icount + 1;    %%%this line means the Jacobian is only calculated once, on the first loop
        end
        
        if NewtonLoopIter>NewtonLoopMaxIter
            disp(strcat('Newton Method interations exceeded MaxIter=',num2str(NewtonLoopMaxIter), ' , loop exit forced.'))
            err=NaN;
        end
    end
    
    
    if isnan(err);
        eta_solved=zeros(N,1);
        disp('KDV solver failed to converge, initial guess set to zero.')
        
    else
        disp('KDV solver converged, initial guess set to KDV solution.')
        if mu <0
            eta_solved=flip(w); % flip soln around in subcritical
        else
            eta_solved=w;
        end
        
    end


    

    
    
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
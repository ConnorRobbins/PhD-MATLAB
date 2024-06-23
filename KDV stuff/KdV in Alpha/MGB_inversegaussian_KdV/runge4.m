    function [x,yout] = runge4(y,dydx,x)
%-- to perform one 4th order Runge-Kutta step.
%     From Numerical Recipes p.706
%     -- One RK step advances the solution from x to (x+h)

    global alpha n h
    
%     -- First step
      for i=1:n
         yt(i)=y(i)+0.5*h*dydx(i);
      end
      [dyt] = eqsystem(x+0.5*h,yt);

%     -- Second step
      for i=1:n
         yt(i)=y(i)+0.5d0*h*dyt(i);
      end
      [dym] = eqsystem(x+0.5*h,yt);

%     -- Third step      
      for i=1:n
        yt=y+h*dym;
        dym=dyt+dym;
      end
%     -- Fourth step
      [dyt] = eqsystem(x+0.5*h,yt);

      for i=1:n
         yout(i)=y(i)+h*(dydx(i)+2*dym(i)+dyt(i))/6;
      end

      x=x+h;

      end

    function [x,yout] = runge4(y,x,h,equationsystem)
%-- to perform one 4th order Runge-Kutta step.
%     EDITED FROM MARK'S CODE:
%     From Numerical Recipes p.706
%     -- One RK step advances the solution from x to (x+h)

    
%     -- First step
      [dydx] = equationsystem(x,y);
      yt=y+0.5*h*dydx;
      [dyt] = equationsystem(x+0.5*h,yt);

%     -- Second step
      yt=y+0.5d0*h*dyt;
      [dym] = equationsystem(x+0.5*h,yt);

%     -- Third step
      yt=y+h*dym;
      dym=dyt+dym;

%     -- Fourth step
      [dyt] = equationsystem(x+0.5*h,yt);


      yout=y+h*(dydx+2*dym+dyt)/6;
      x=x+h;

      end

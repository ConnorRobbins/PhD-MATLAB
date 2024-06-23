function [dydx] = eqsystem(x,y)
     
    dydx(1) = y(2);
    dydx(2) =-sin(y(1));
    
end
function [dydx] = eqsystem2(x,y)
     
    global alpha 
              
    ff = 2*(2*x^2-1)*exp(-x^2) + exp(-2*x^2);
    
    dydx(1) = y(2);
    dydx(2) = alpha*ff - y(1)^2;
    
end
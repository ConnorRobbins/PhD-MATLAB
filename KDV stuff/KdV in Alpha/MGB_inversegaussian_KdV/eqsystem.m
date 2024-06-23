function [F] = eqsystem(x,y)
     
    global alpha n h
              
    ff = 2*(2*x^2-1)*exp(-x^2) + exp(-2*x^2);
    
    F(1) = y(2);
    F(2) = alpha*ff - y(1)^2;
    
end
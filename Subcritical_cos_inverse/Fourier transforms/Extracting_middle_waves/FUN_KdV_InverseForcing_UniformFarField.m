function [forcing] = FUN_KdV_InverseForcing_UniformFarField(eta,x,Froude)
%KDV_INVERSEFORCING_UNIFORMFARFIELD Calculates the inverse forcing with a
%uniform far field condition from the fKdV 

del_x=x(2)-x(1);
eta_xx=[0,eta(3:end)-2*eta(2:end-1)+eta(1:end-2),0]/del_x^2;
forcing=2*(Froude-1)*eta - (eta_xx/3)  - 1.5*eta.^2;

end


clc
%clear

L_width = 15; %Phi range on either side of the origin
N_points = 21;
Froude = 1.1;
iload = 0;
isave = 0;
ijac = 1;
Pressure = 0;
Topography = 1;

amplitude=0.002;
b_width=sqrt(1);   %sqrt(ben's amplitude)


N_s_points=N_points;
N_t_points=N_points;


main_forward;

Forward_GN;
jac_GN=full(jac_GN);


%% Trying to match the JACs
jac_GN2=jac_GN(1:end-1,2:end);

jac_diff=jac_GN2-jac_N;
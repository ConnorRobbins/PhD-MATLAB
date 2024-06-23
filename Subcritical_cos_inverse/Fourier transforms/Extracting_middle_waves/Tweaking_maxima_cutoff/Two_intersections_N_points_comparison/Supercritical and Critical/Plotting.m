clear

%load('p4_to_1_crit.mat')
load('datacleaned.mat')
%load('Apoint05_N1624_k_half_to_2.mat')
%load('Apoint01_N1624_k_half_to_2.mat')
%load('Apoint05_60_ks_N1624_k_half_to_1.mat')
%load('joined_data.mat')
%load('subcrit_2inter_0p7_to_0p9875.mat')
load('sub_p2_to_1p1.mat')

for i=1:ki
    kdv_freq=kdv_freqs{i};
    kdv_y=kdv_yshift{i};
    nonlinear_freq=nonlinear_freqs{i};
    nonlinear_y=nonlinear_yshift{i};

    k=ks(i);
    
    % Extract c1 mode
    [~,kdv_c1_index]=min(abs(kdv_freq-k));
    [~,nonlinear_c1_index]=min(abs(nonlinear_freq-k));
    kdv_c1_power(i)=kdv_y(kdv_c1_index);
    nonlinear_c1_power(i)=nonlinear_y(nonlinear_c1_index);


     % Extract c2 mode
    [~,kdv_c2_index]=min(abs(kdv_freq-2*k));
    [~,nonlinear_c2_index]=min(abs(nonlinear_freq-2*k));
    kdv_c2_power(i)=kdv_y(kdv_c2_index);
    nonlinear_c2_power(i)=nonlinear_y(nonlinear_c2_index);


    [~,nonlinear_c3_index]=min(abs(nonlinear_freq-3*k));
    nonlinear_c3_power(i)=nonlinear_y(nonlinear_c3_index);

end


figure(61); clf; hold on;
plot(ks,nonlinear_c1_power,'-xk',DisplayName='non c1')
plot(ks,nonlinear_c2_power,'--xk',DisplayName='non c2')
% plot(ks,kdv_c1_power,'-xr')
% plot(ks,kdv_c2_power,'--xr')


plot(ks,ones(size(ks))*3*(amplitude^2)/4,'--b',DisplayName='kdv c2')
plot(ks,abs(amplitude*(2*(Froude-1) +(ks.^2)/3)),'-xb',DisplayName='kdv c1')
legend

%legend('nonlin c1','nonlin c2','numeric kdv c1','numeric kdv c2','analytic kdv c1','analytic kdv c2')
%plot(ks,nonlinear_c3_power,':m')





figure(62); clf; hold on;
plot(ks,log10(nonlinear_c1_power),'-k')
plot(ks,log10(nonlinear_c2_power),'-r')
%plot(ks,log10(nonlinear_c3_power),'-c')
legend('log(c1)','log(c2)')


figure(63); clf; hold on; 
[straight_c1,m1]=FUN_simple_straight_line_fit(ks,log10(nonlinear_c1_power));
[straight_c2,m2]=FUN_simple_straight_line_fit(ks,log10(nonlinear_c2_power));
plot(ks,log10(nonlinear_c1_power),'-k')
plot(ks,straight_c1,'--g')
plot(ks,log10(nonlinear_c2_power),'-r')
plot(ks,straight_c2,'--c')
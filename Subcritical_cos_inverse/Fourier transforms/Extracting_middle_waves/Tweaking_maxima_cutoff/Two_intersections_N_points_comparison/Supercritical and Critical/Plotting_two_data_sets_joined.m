clear


%file1='Apoint05_k_point5_to_1point94.mat';
file1='datacleaned.mat';
%file1='subcrit_2inter_0p7_to_0p9925';
%file1='Apoint05_N1624_k_half_to_2.mat';


%file2="p6_to_2";
%file2='Apoint05_60_ks_N1624_k_1_to_2.mat';
%file2="datacleaned";
%file2='Apoint05_60_ks_N1624_k_half_to_1.mat';
%file2="Apoint05_N1624_k_half_to_2.mat";
file2='subcrit_2inter_0p7_to_1';




%% Do it for first data set
load(file1)




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
plot(ks,nonlinear_c1_power,'-k')
plot(ks,nonlinear_c2_power,'--k')
% plot(ks,kdv_c1_power,'-xr')
% plot(ks,kdv_c2_power,'--xr')


% plot(ks,ones(size(ks))*3*(amplitude^2)/4,'--b')
% plot(ks,abs(amplitude*(2*(Froude-1) +(ks.^2)/3)),'-xb')

legend('nonlin c1','nonlin c2','numeric kdv c1','numeric kdv c2','analytic kdv c1','analytic kdv c2')
%plot(ks,nonlinear_c3_power,':m')





figure(62); clf; hold on;
plot(ks,log10(nonlinear_c1_power),'-k')
plot(ks,log10(nonlinear_c2_power),'-xr')
legend('log(c1)','log(c2)')




%% Do it for second data set
clearvars -except file2
load(file2)




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


figure(61); 
plot(ks,nonlinear_c1_power,'-xk')
plot(ks,nonlinear_c2_power,'--xk')
% plot(ks,kdv_c1_power,'-xr')
% plot(ks,kdv_c2_power,'--xr')


% plot(ks,ones(size(ks))*3*(amplitude^2)/4,'--b')
% plot(ks,abs(amplitude*(2*(Froude-1) +(ks.^2)/3)),'-xb')

legend('nonlin c1','nonlin c2','numeric kdv c1','numeric kdv c2','analytic kdv c1','analytic kdv c2')
%plot(ks,nonlinear_c3_power,':m')





figure(62); 
plot(ks,log10(nonlinear_c1_power),'-xk')
plot(ks,log10(nonlinear_c2_power),'-xr')
legend('log(c1)','log(c2)')

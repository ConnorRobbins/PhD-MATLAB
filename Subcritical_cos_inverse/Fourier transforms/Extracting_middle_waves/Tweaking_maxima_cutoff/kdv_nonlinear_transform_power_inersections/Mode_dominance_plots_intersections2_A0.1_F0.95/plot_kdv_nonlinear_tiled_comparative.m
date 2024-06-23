clear

s=load('right_intersection.mat');
right_freq_kdv=s.pos_frequencies_kdv_middle_only;
right_freq_nonlinear=s.pos_frequencies_nonlinear_middle_only;
right_yshift_kdv=s.pos_yshift_kdv_middle_only;
right_yshift_nonlinear=s.pos_yshift_nonlinear_middle_only;
right_wave_k=s.wave_k;

s=load('c1_equals_0.mat');
central_freq_kdv=s.pos_frequencies_kdv_middle_only;
central_freq_nonlinear=s.pos_frequencies_nonlinear_middle_only;
central_yshift_kdv=s.pos_yshift_kdv_middle_only;
central_yshift_nonlinear=s.pos_yshift_nonlinear_middle_only;
central_wave_k=s.wave_k;

s=load('left_intersection.mat');
left_freq_kdv=s.pos_frequencies_kdv_middle_only;
left_freq_nonlinear=s.pos_frequencies_nonlinear_middle_only;
left_yshift_kdv=s.pos_yshift_kdv_middle_only;
left_yshift_nonlinear=s.pos_yshift_nonlinear_middle_only;
left_wave_k=s.wave_k;

%% 
xticklabs={'$0$','$k$','$2k$','$3k$','$4k$'};
xleft_ticks=[0,left_wave_k,2*left_wave_k,3*left_wave_k,4*left_wave_k];
xcentral_ticks=[0,central_wave_k,2*central_wave_k,3*central_wave_k,4*central_wave_k];
xright_ticks=[0,right_wave_k,2*right_wave_k,3*right_wave_k,4*right_wave_k];


xmulti=0.4;
ylims=([0,0.016]);
circlesize=4;

%set(groot,'defaultAxesTickLabelInterpreter','latex'); 


figure(15); clf; hold on;
fi=tiledlayout(1,3);
nexttile
stem(left_freq_kdv,left_yshift_kdv,MarkerSize=circlesize,Color='k')
hold on
stem(left_freq_nonlinear,left_yshift_nonlinear,MarkerSize=circlesize,Color='r')
xlim([-xmulti*left_wave_k,(4+xmulti)*left_wave_k])
ylim(ylims)
xticks(xleft_ticks)
xticklabels(xticklabs)


nexttile
stem(central_freq_kdv,central_yshift_kdv,MarkerSize=circlesize,Color='k')
hold on
stem(central_freq_nonlinear,central_yshift_nonlinear,MarkerSize=circlesize,Color='r')
xlim([-xmulti*central_wave_k,(4+xmulti)*central_wave_k])
ylim(ylims)
xticks(xcentral_ticks)
xticklabels(xticklabs)



nexttile
stem(right_freq_kdv,right_yshift_kdv,MarkerSize=circlesize,Color='k')
hold on
stem(right_freq_nonlinear,right_yshift_nonlinear,MarkerSize=circlesize,Color='r')
xlim([-xmulti*right_wave_k,(4+xmulti)*right_wave_k])
ylim(ylims)
xticks(xright_ticks)
xticklabels(xticklabs)


xlabel(fi,'Frequency')
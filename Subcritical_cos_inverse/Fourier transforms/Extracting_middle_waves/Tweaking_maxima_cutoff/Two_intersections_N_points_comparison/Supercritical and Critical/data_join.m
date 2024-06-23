clear
file2=load('subcrit_2inter_0p7_to_1.mat');
file1=load('datacleaned.mat');
kdv_freqs={file1.kdv_freqs{1,:},file2.kdv_freqs{1,:}};
kdv_yshift={file1.kdv_yshift{1,:},file2.kdv_yshift{1,:}};
nonlinear_freqs={file1.nonlinear_freqs{1,:},file2.nonlinear_freqs{1,:}};
nonlinear_yshift={file1.nonlinear_yshift{1,:},file2.nonlinear_yshift{1,:}};
ranks=[file1.ranks,file2.ranks];
Ys_store={file1.Ys_store{1,:},file2.Ys_store{1,:}};
Yt_kdv_store={file1.Yt_kdv_store{1,:},file2.Yt_kdv_store{1,:}};
Yt_nonlinear_store={file1.Yt_nonlinear_store{1,:},file2.Yt_nonlinear_store{1,:}};

ks=[file1.ks,file2.ks];
ki=file1.ki+file2.ki;
Froude=file1.Froude;
amplitude=file1.amplitude;
clear file1 file2
save('joined_data.mat')







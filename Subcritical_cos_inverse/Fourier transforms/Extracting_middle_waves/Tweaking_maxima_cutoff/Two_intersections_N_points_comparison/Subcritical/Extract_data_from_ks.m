clear

ki=200;


ks=zeros(1,ki);
ranks=zeros(1,ki);



for i=1:ki
    loadstring=strcat(num2str(ki),'K',num2str(i));
    s=load(loadstring);
    ks(i)=s.wave_k;
    kdv_freqs{i}=s.pos_frequencies_kdv_middle_only;
    nonlinear_freqs{i}=s.pos_frequencies_nonlinear_middle_only;
    kdv_yshift{i}=s.pos_yshift_kdv_middle_only;
    nonlinear_yshift{i}=s.pos_yshift_nonlinear_middle_only;
    Ys_store{i}=s.Ys;
    Yt_kdv_store{i}=s.Yt_kdv;
    Yt_nonlinear_store{i}=s.Yt_nonlinear;
    ranks(i)=s.rankSVD;
end

amplitude=s.amplitude;
Froude=s.Froude;

clear s loadstring
save('datacleaned.mat')


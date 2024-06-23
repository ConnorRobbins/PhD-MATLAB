Ns=[81,141,241,341,441,841,1241,1641];
SVDranks=[81,141,119,116,111,104,99,91];


NNs=numel(Ns);

filename="unforced_solitary_F1p1_N";
maxvalstore=zeros(1,NNs);




for i=1:NNs
    load(filename+num2str(Ns(i)));
    [Phi_t,Phi_s,Theta_inverse,~,~,~,Theta_b_matrix_inverse,~,~,~,~,~] = svd_Step_fun(L,Froude,N,N,P,Ys_Newton);
    Theta_bottom_inverse=Theta_b_matrix_inverse(:,SVDranks(i));
    [Yt_inverse,Ys_inverse,Xt_inverse,Xs_inverse] = variables_after_SVD(L,N,N,Theta_bottom_inverse,Theta_inverse);
    maxvalstore(i)=max(abs(Yt_inverse));
end



figure(4678); clf; hold on;
plot(Ns,maxvalstore,'-xr')
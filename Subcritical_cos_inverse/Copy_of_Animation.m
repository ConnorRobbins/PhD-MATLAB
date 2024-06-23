%%% data must have same number of steps
data_1=load('looping_topo_centre_5_25_200steps_rankSVD400_amp0.05.mat');
data_2=load('looping_topo_centre_5_25_200steps_rankSVD400_amp0.1.mat');

topo_centres=data_1.topo_centres;

Ys_store_1=data_1.Ys_store;
Yt_store_1=data_1.Yt_store;
Phi_s_1=data_1.Phi_s;
Phi_t_1=data_1.Phi_t;

Ys_store_2=data_2.Ys_store;
Yt_store_2=data_2.Yt_store;
Phi_s_2=data_2.Phi_s;
Phi_t_2=data_2.Phi_t;

total_frames=length(topo_centres);
Ys_max=max(max(max(Ys_store_1)),max(max(Ys_store_2)));
Ys_min=min(min(min(Ys_store_1)),min(min(Ys_store_2)));
Yt_max=max(max(max(Yt_store_1)),max(max(Yt_store_2)));
Yt_min=min(min(min(Yt_store_1)),min(min(Yt_store_2)));




%% Video attempt

figure(1)
clf
subplot(2,1,1)
ylabel('Ys')
subplot(2,1,2)
ylabel('Yt')

% xlabel('$\phi$','Interpreter','Latex','FontSize',16)
% ylabel('$\theta_b$','Interpreter','Latex','FontSize',18,'Rotation',0)



for k=1:total_frames
    subplot(2,1,1)
    plot(Phi_s_1,Ys_store_1(:,k),'-r',Phi_s_2,Ys_store_2(:,k),'-b')
    ylim([Ys_min,Ys_max])
    subplot(2,1,2)
    plot(Phi_t_1,Yt_store_1(:,k),'-r',Phi_t_2,Yt_store_2(:,k),'-b')
    ylim([Yt_min,Yt_max])
    Frames(k)=getframe(gcf);
    drawnow
end



% create the video writer with 1 fps
writerObj = VideoWriter('myVideo_two_sets.avi');
writerObj.FrameRate = 10;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for k=1:total_frames
    % convert the image to a frame
    frame = Frames(k) ;
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);

implay('myVideo_two_sets.avi')
%load('looping_topo_centre.mat')
load('looping_topo_centre_5_25_200steps_rankSVD400_amp0.1.mat')

total_frames=length(topo_centres);
Ys_max=max(max(Ys_store));
Ys_min=min(min(Ys_store));
Yt_max=max(max(Yt_store));
Yt_min=min(min(Yt_store));



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
    plot(Phi_s,Ys_store(:,k))
    ylim([Ys_min,Ys_max])
    subplot(2,1,2)
    plot(Phi_s,Yt_store(:,k))
    ylim([Yt_min,Yt_max])
    Frames(k)=getframe(gcf);
    drawnow
end



% create the video writer with 1 fps
writerObj = VideoWriter('myVideo.avi');
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

implay('myVideo.avi')
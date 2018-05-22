function B_mode_disp(x,t,fr)
%x is 3dimensional matrix containing data
%t is a number or vector containing the time to view the B mode image

if length(t) == 1
t1 = x(:,:,t);
t2 = squeeze(t1);
imagesc(t2)
else
%     vidObj = VideoWriter('bmode.avi');
%     vidObj.FrameRate=fr;
%     open(vidObj);
    
for i = t
    t1 = x(:,:,i);
    t2 = squeeze(t1);
    imagesc(t2)
    colormap(jet)
    f(i-min(t)+1) = getframe;
 %   writeVideo(vidObj,f(i))
end

%movie(f,1,fr)
% implay(vidObj);
%close(vidObj);
end

    
    
   
  
    

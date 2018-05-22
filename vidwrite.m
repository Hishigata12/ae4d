function vidwrite(y,AxZ,AxX,z,x,t,c)
v = VideoWriter('test.avi');
open(v);
y = squeeze(y);
if exist('x') == 0
    x(1) = 1;
    x(2) = size(y,1);
end
if exist('z') == 0
    z(1) = 1;
    z(2) = size(y,2);
end

if exist('t') == 0
    t(1) = 1;
    t(2) = size(y,3); 
end
figure;
%y2 = real(20*log10(y/max(max(max(y)))));
for i = t(1):t(2)
    imagesc(AxX(x(1):x(2)),AxZ(z(1):z(2)),medfilt2(y(x(1):x(2),z(1):z(2),i)',[5 5]));
    title(['t = ' num2str(i)])
    colormap('hot');
    if exist('c','var')
        caxis(c)
    end
    frame = getframe;
    writeVideo(v,frame);
end
close
close(v)
end
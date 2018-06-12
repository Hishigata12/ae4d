function vidwrite(param,ax,Xfilt,handles)
file = handles.savefigname.String;
path = handles.savefolder.String;
% v = VideoWriter([path file]);
v = VideoWriter(file);
v.FrameRate = str2double(handles.framerate.String);
open(v);


xR = str2num(handles.xR.String);
if length(xR) == 1
    xR = [xR xR];
end
yR = str2num(handles.yR.String);
if length(yR) == 1
    yR = [yR yR];
end
zR = str2num(handles.zR.String);
if length(zR) == 1
    zR = [zR zR];
end
tR = str2num(handles.tR.String);
if length(tR) == 1
    tR = [tR tR];
end
aeR = str2num(handles.aeR.String);
dims = size(Xfilt);

q.x = 1:dims(1);
q.y = 1:dims(2);
q.z = 1:dims(3);
if length(dims) == 4
    q.t = 1:dims(4);
else
    q.t = 1;
end
xInd = q.x(find(ax.x >= xR(1)):find(ax.x >= xR(2)));
yInd = q.y(find(ax.y >= yR(1)):find(ax.y >= yR(2)));
zInd = q.z(find(ax.depth >= zR(1)):find(ax.depth >= zR(2)));
tInd = q.t(find(ax.stime >= tR(1)):find(ax.stime >= tR(2)));

Xfilt = squeeze(Xfilt(xInd,yInd,zInd,tInd));

if handles.plotbox2.Value == 1
    n = length(tInd);
    p = 't';
     Ind = tInd;
elseif handles.plotbox2.Value == 2
    n = length(tInd);
    p = 't';
     Ind = tInd;
elseif handles.plotbox2.Value == 3
    n = length(zInd);
    p = 'z';
    Ind = zInd;
end

if handles.hotcold.Value == 1
    h = hotcoldDB;
elseif handles.graybox.Value == 1
    h = 'gray';
else
    h = 'hot';
end

for i = 1:n
    if handles.use_ext_fig.Value == 1
        imshow(squeeze(Xfilt(:,:,i))')
    else
    imagesc(ax.x(xInd),ax.depth(zInd),squeeze(Xfilt(:,:,i))')
    end
    title([p ' = ' num2str(Ind(i))])   
    colormap(h);
    if ~isempty(aeR)
        caxis(aeR)
    end
    frame = getframe;
    writeVideo(v,frame);
end

close(v)
% 
% if exist('x') == 0
%     x(1) = 1;
%     x(2) = size(y,1);
% end
% if exist('z') == 0
%     z(1) = 1;
%     z(2) = size(y,2);
% end
% 
% if exist('t') == 0
%     t(1) = 1;
%     t(2) = size(y,3); 
% end
% figure;
% %y2 = real(20*log10(y/max(max(max(y)))));
% for i = t(1):t(2)
%     imagesc(AxX(x(1):x(2)),AxZ(z(1):z(2)),medfilt2(y(x(1):x(2),z(1):z(2),i)',[5 5]));
%     title(['t = ' num2str(i)])
%     colormap('hot');
%     if exist('c','var')
%         caxis(c)
%     end
%     frame = getframe;
%     writeVideo(v,frame);
% end
% close
% close(v)
% end
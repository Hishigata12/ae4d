function vidwrite(Xfilt,handles,aeR,h)
file = handles.savefigname.String;
path = handles.savefolder.String;
% v = VideoWriter([path file]);
v = VideoWriter(file);
v.FrameRate = str2double(handles.framerate.String);
open(v);

% 
% xR = str2num(handles.xR.String);
% if length(xR) == 1
%     xR = [xR xR];
% end
% yR = str2num(handles.yR.String);
% if length(yR) == 1
%     yR = [yR yR];
% end
% zR = str2num(handles.zR.String);
% if length(zR) == 1
%     zR = [zR zR];
% end
% tR = str2num(handles.tR.String);
% if length(tR) == 1
%     tR = [tR tR];
% end
% if exist('handles.aeR','var')
% aeR = str2num(handles.aeR.String);
% elseif exist('handles.cmin')
%     aeR = [str2num(handles.cmin.String) str2num(handles.cmax.String)];
% end
% dims = size(Xfilt);
% 
% q.x = 1:dims(1);
% q.y = 1:dims(2);
% q.z = 1:dims(3);
% if length(dims) == 4
%     q.t = 1:dims(4);
% else
%     q.t = 1;
% end
% xInd = q.x(find(ax.x >= xR(1)):find(ax.x >= xR(2)));
% yInd = q.y(find(ax.y >= yR(1)):find(ax.y >= yR(2)));
% zInd = q.z(find(ax.depth >= zR(1)):find(ax.depth >= zR(2)));
% tInd = q.t(find(ax.stime >= tR(1)):find(ax.stime >= tR(2)));
%   xP = str2double(handles.xP.String);
%     yP = str2double(handles.yP.String);
%     zP = str2double(handles.zP.String);
%     tP = str2double(handles.tP.String);
%     xpoint = find(ax.x >= xP,1);
%     ypoint = find(ax.y >= yP,1);
%     zpoint = find(ax.depth >= zP,1);
%     tpoint = find(ax.stime >= tP,1);


% Xfilt = squeeze(Xfilt(xInd,yInd,zInd,tInd));
n = size(Xfilt,3);
p = 'n';


% if handles.plotbox2.Value == 1
%     n = length(tInd);
%     p = 't';
%      Ind = tInd;
%      Xfilt = squeeze(Xfilt(xInd,ypoint,zInd,tInd));
% elseif handles.plotbox2.Value == 2
%     n = length(tInd);
%     p = 't';
%      Ind = tInd;
%      Xfilt = squeeze(Xfilt(xpoint,yInd,zInd,tInd));
% elseif handles.plotbox2.Value == 3
%     n = length(zInd);
%     p = 'z';
%     Ind = zInd;
%     Xfilt = squeeze(Xfilt(xInd,yInd,zInd,tpoint));
% elseif handles.plotbox2.Value == 4
%       n = length(tInd);
%     p = 't';
%      Ind = tInd;
%      Xfilt = squeeze(Xfilt(xInd,yInd,zpoint,tInd));
% end


% if handles.hotcold.Value == 1
%     if handles.bbdb.Value == 1
%         h = hotcoldDB;
%     else
%         h = hotcold;
%     end
% elseif handles.graybox.Value == 1
%     h = 'gray';
% else
%     h = 'hot';
% end
if handles.showlf.Value
    LF = evalin('base','LF');
    chan = str2double(handles.LF_chan.String);
    param = evalin('base','param');
    lf_ax = linspace(0,param.Duration,length(LF));
    LF = LF(:,chan);
    
    %build time axis stuff
    if handles.use_chop.Value
        ax = evalin('base','ax_c');
    else
        ax = evalin('base','ax');
    end
    tR = str2num(handles.tR.String);
    if length(tR) == 1
        tR = [tR tR];
    end
    
    dims = size(Xfilt);
    
    if length(dims) == 3
        q.t = 1:dims(3);
    else
        q.t = 1;
    end
    
    tInd = q.t(find(ax.stime >= tR(1)):find(ax.stime >= tR(2)));
    
    tP = str2double(handles.tP.String);
    
    tpoint = find(ax.stime >= tP,1);
    
    lfInd = find(lf_ax >= tR(1)):find(lf_ax >= tR(2));
    lfdif = length(lfInd)/length(tInd);
    
     for k = tInd
     plot(lf_ax(lfInd),LF(lfInd),'w')
            hold('on')
%             xlim([tR(1) tR(2)]);
            %  plot(handles.axes4,lf_ax(lfInd(round(k*lfdif))),LF(lfInd(round(k*lfdif))),'ro','MarkerFaceColor','r')
%             plot(handles.axes4,lf_ax(round(k*lfdif)),LF(round(k*lfdif)),'ro','MarkerFaceColor','r')
            plot(lf_ax(lfInd(round(k*lfdif))),LF(lfInd(round(k*lfdif))),'ro','MarkerFaceColor','r')
            hold('off');
            set(gca,'Color','k');
             frame = getframe;
    writeVideo(v,frame);
     end
else
    for i = 1:n
        %     if handles.use_ext_fig.Value == 1
        %     if i == 1
        %         g = axes;
        imshow(Xfilt(:,:,i)');
        %     else
        %         figure(145)
        %         squeeze(Xfilt(:,:,i))';
        %         %  pic.CData = squeeze(Xfilt(:,:,i))';
        %     end
        %     else
        %     imagesc(squeeze(Xfilt(:,:,i))')
        %     end
        title([p ' = ' num2str(i)])
        colormap(h);
        if ~isempty(handles.movietext.String)
            text(1,15,handles.movietext.String,'Color','white','FontSize',14)
        end
        if ~isempty(aeR)
            caxis(aeR)
        end
        frame = getframe;
        writeVideo(v,frame);
    end
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